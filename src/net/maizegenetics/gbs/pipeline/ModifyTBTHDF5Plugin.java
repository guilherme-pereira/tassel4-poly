/*
 * ModifyTBTHDF5Plugin
 */
package net.maizegenetics.gbs.pipeline;

import cern.colt.list.IntArrayList;

import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteHDF5TaxaGroups;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteHDF5TagGroups;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

import java.awt.Frame;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

/**
 *  This pipeline modifies TagsByTaxa HDF5 file with data organized by taxa.
 * It can:
 * 1.  Create an empty TBT.
 * 2.  Merge two TBT
 * 3.  Combined similarly named taxa
 * 4.  Pivot a taxa TBT to a tag TBT
 *
 * @author ed 
 */
public class ModifyTBTHDF5Plugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ModifyTBTHDF5Plugin.class);
    private ArgsEngine myArgsEngine = null;
    private String myAdditionalTBTHDF5 = null;
    private String myTargetTBTHDF5 = null;
    private String myOutputTransposeTagTBTHDF5 = null;
//    private String myOutputLogFile = null;
    private boolean combineTaxa = false;
    private Tags myMasterTags = null;
    private static int maxGoodReads = 500000000; // maximum number of good barcoded reads expected in a fastq file
    IntArrayList[] taxaReads;
    int[] readsPerSample, mappedReadsPerSample;
    int goodBarcodedReads = 0, allReads = 0, goodMatched = 0;
    HashMap<String, Integer> taxaNameToIndices;
    TagsByTaxaByteHDF5TaxaGroups targetTBT = null;

    public ModifyTBTHDF5Plugin() {
        super(null, false);
    }

    public ModifyTBTHDF5Plugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        if (myMasterTags != null) {
            targetTBT = new TagsByTaxaByteHDF5TaxaGroups(myMasterTags, myTargetTBTHDF5);
        } else {
            targetTBT = new TagsByTaxaByteHDF5TaxaGroups(myTargetTBTHDF5);
        }
        if (myAdditionalTBTHDF5 != null) {
            addAllTaxaToNewHDF5(myAdditionalTBTHDF5);
        }
        if (combineTaxa) {
            combineTaxaHDF5();
        }
        if (myOutputTransposeTagTBTHDF5 != null) {
            TagsByTaxaByteHDF5TagGroups tranTBT = new TagsByTaxaByteHDF5TagGroups(targetTBT, myOutputTransposeTagTBTHDF5);
        }
        targetTBT.getFileReadyForClosing();
        targetTBT = null;
        System.gc();
        return null;
    }

    private void printUsage() {
        myLogger.info(
                "\n\n\nThe ModifyTBTHDF5Plugin accepts the following arguments:\n"
                + "-o  Target TBT HDF5 (*tbt.h5) file to be modified\n"
                + "Depending on the modification that you wish to make, one of either:\n"
                + "    -i  TBT HDF5 (*tbt.h5) file containing additional taxa to be added to the target TBT HDF5 file\n"
                + "    -c  Merge taxa in the target TBT HDF5 file with same LibraryPrepID\n"
                + "    -p  Pivot (transpose) the target TBT HDF5 file into a tag-optimized orientation\n"
//                + "For creating an emptry TBT HDF4, one of either:\n"
//                + "    -t  Tag count file, OR A\n"
//                + "    -m  Physical map file containing alignments\n"
//                + "-L  Output log file \n"
                +"\n\n\n");
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-i", "--addition-file", true);
            myArgsEngine.add("-o", "--target-HDF5", true);
//            myArgsEngine.add("-L", "--outputlogfile", true);
            myArgsEngine.add("-c", "--combine-taxa", false);
            myArgsEngine.add("-p", "--taghdf-file", true);
//            myArgsEngine.add("-t", "--tagcount-file", true);
//            myArgsEngine.add("-m", "--physical-map", true);
        }
        myArgsEngine.parse(args);
        if (myArgsEngine.getBoolean("-o")) {
            myTargetTBTHDF5 = myArgsEngine.getString("-o");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an target file (option -o).");
        }
        if (myArgsEngine.getBoolean("-i")) {
            myAdditionalTBTHDF5 = myArgsEngine.getString("-i");
        }
        if (myArgsEngine.getBoolean("-p")) {
            myOutputTransposeTagTBTHDF5 = myArgsEngine.getString("-p");
        }
        if (myArgsEngine.getBoolean("-c")) {
            combineTaxa = true;
        }
//        if (myArgsEngine.getBoolean("-L")) {
//            myOutputLogFile = myArgsEngine.getString("-L");
//        } else {
//            printUsage();
//            throw new IllegalArgumentException("Please specify a log file (option -L).");
//        }
//        // Create Tags object from tag count file with option -t, or from TOPM file with option -m
//        if (myArgsEngine.getBoolean("-t")) {
//            if (myArgsEngine.getBoolean("-m")) {
//                printUsage();
//                throw new IllegalArgumentException("Options -t and -m are mutually exclusive.");
//            }
//            myMasterTags = new TagCounts(myArgsEngine.getString("-t"), FilePacking.Bit);
//        } else if (myArgsEngine.getBoolean("-m")) {
//            if (myArgsEngine.getBoolean("-t")) {
//                printUsage();
//                throw new IllegalArgumentException("Options -t and -m are mutually exclusive.");
//            }
//            myMasterTags = new TagsOnPhysicalMap(myArgsEngine.getString("-m"), true);
//        }
    }

    private boolean addAllTaxaToNewHDF5(String addTBTName) {
        TagsByTaxaByteHDF5TaxaGroups addTBT = new TagsByTaxaByteHDF5TaxaGroups(addTBTName);
        for (int i = 0; i < addTBT.getTagCount(); i++) {
            if (!Arrays.equals(targetTBT.getTag(i), addTBT.getTag(i))) {
                System.err.println("Tags are not the same for the two TBT file.  They cannot be combined.");
                System.exit(0);
            }
        }
        for (int i = 0; i < addTBT.getTaxaCount(); i++) {
            String name = addTBT.getTaxaName(i);
            byte[] states = addTBT.getReadCountDistributionForTaxon(i);
            targetTBT.addTaxon(name, states);
        }
        addTBT.getFileReadyForClosing();
        return true;
    }

    private boolean combineTaxaHDF5() {
        TreeMap<String, ArrayList<String>> combineTaxa = new TreeMap<String, ArrayList<String>>();
        for (String tn : targetTBT.getTaxaNames()) {
            String ptn = parseTaxaName(tn, "#");
            ArrayList<String> taxaList = combineTaxa.get(ptn);
            if (taxaList == null) {
                combineTaxa.put(ptn, taxaList = new ArrayList<String>());
            }
            taxaList.add(tn);
        }
        for (ArrayList<String> taxaList : combineTaxa.values()) {
            if (taxaList.size() < 2) {
                continue;
            }
            byte[] di = new byte[targetTBT.getTagCount()];
            String ptn = parseTaxaName(taxaList.get(0), "" + taxaList.size());
            for (int i = 0; i < taxaList.size(); i++) {
                int j = targetTBT.getIndexOfTaxaName(taxaList.get(i));
                byte[] dj = targetTBT.getReadCountDistributionForTaxon(j);
                for (int k = 0; k < dj.length; k++) {
                    int ts = di[k] + dj[k];
                    if (ts > Byte.MAX_VALUE) {
                        di[k] = Byte.MAX_VALUE;
                        // System.out.println(di[k]+"+"+dj[k]);
                    } else {
                        di[k] = (byte) ts;
                    }
                }
            }
            targetTBT.addTaxon(ptn, di);
            for (String tn : taxaList) {
                targetTBT.deleteTaxon(tn);
            }
        }
        return true;
    }

    private String parseTaxaName(String tn, String cnt) {
        String[] s = tn.split(":");
        return s[0] + ":MRG:" + cnt + ":" + s[3];
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

    public static void main(String[] args) {
        String inTBTFile = "/Users/edbuckler/SolexaAnal/GBS/build20120110/tbt/434GFAAXX_s_3.tbt.byte";
        TagsByTaxa inTBT = new TagsByTaxaByte(inTBTFile, FilePacking.Byte);
//        TagsByTaxaByteHDF5TaxaGroups tHDF5 = new TagsByTaxaByteHDF5TaxaGroups(inTBT, "/Users/edbuckler/SolexaAnal/GBS/test/433s3.h5");
//        Tags theTags=new TagCountMutable(inTBT, inTBT.getTagCount());
//        TagsByTaxaByteHDF5TaxaGroups rHDF5 = new TagsByTaxaByteHDF5TaxaGroups(theTags, "/Users/edbuckler/SolexaAnal/GBS/test/r433s3.h5");
//        for (int i=0; i<tHDF5.getTaxaCount(); i++) {
//            byte[] d=tHDF5.getReadCountDistributionForTaxon(i);
//            String newName=tHDF5.getTaxaName(i).replace(":3:", ":4:");
//            rHDF5.addTaxon(newName, d);
//        }
//        tHDF5.getFileReadyForClosing();
//        rHDF5.getFileReadyForClosing();
//        
//        args = new String[] {
//            "-i", "/Users/edbuckler/SolexaAnal/GBS/test/r433s3.h5",
//            "-o", "/Users/edbuckler/SolexaAnal/GBS/test/433s3.h5",
//            "-L", "/Users/edbuckler/SolexaAnal/GBS/test/testFQSeqTBT.log",
//            //"-c",
//            "-p","/Users/edbuckler/SolexaAnal/GBS/test/t433s3.h5",
//        };
//
//        ModifyTBTHDF5Plugin testClass = new ModifyTBTHDF5Plugin();
//        testClass.setParameters(args);
//        testClass.performFunction(null);
//        
//         rHDF5 = new TagsByTaxaByteHDF5TaxaGroups("/Users/edbuckler/SolexaAnal/GBS/test/433s3.h5");
//         tHDF5 = new TagsByTaxaByteHDF5TaxaGroups("/Users/edbuckler/SolexaAnal/GBS/test/r433s3.h5");
//
//        for (int i=0; i<rHDF5.getTaxaCount(); i++) {
//            byte[] di=rHDF5.getReadCountDistributionForTaxon(i);
//            int sumi=0;
//            for (byte b : di) {sumi+=b;}
//            System.out.println(rHDF5.getTaxaName(i)+":"+sumi);
//            for (int j = 0; j < tHDF5.getTaxaCount(); j++) {
//                byte[] dj=tHDF5.getReadCountDistributionForTaxon(j);
//                if(Arrays.equals(di, dj)) {
//                    int sumj=0;
//                    for (byte b : dj) {sumj+=b;}
//                    System.out.println(rHDF5.getTaxaName(i)+":"+sumi+"="+tHDF5.getTaxaName(j)+":"+sumj);
//                }
//            }
//        }
        TagsByTaxaByteHDF5TagGroups transHDF5 = new TagsByTaxaByteHDF5TagGroups("/Users/edbuckler/SolexaAnal/GBS/test/t433s3.h5");
        int same = 0, diff = 0, count = 0;
        long time = System.currentTimeMillis();
        int tags = 9;
        for (int i = 0; i < 1000000; i += 11) {
            int taxon = i % inTBT.getTaxaCount();
            // taxon=15;
            if (i % 550 == 0) {
                tags = i % inTBT.getTagCount();
            }

            int newTaxonIndex = transHDF5.getIndexOfTaxaName(inTBT.getTaxaName(taxon));
            if (inTBT.getReadCountForTagTaxon(tags, taxon) == transHDF5.getReadCountForTagTaxon(tags, newTaxonIndex)) {
                same++;
            } else {
                diff++;
            }
            count++;
        }
        System.out.printf("Same %d Diff %d %n", same, diff);
        long duration = System.currentTimeMillis() - time;
        double rate = (double) duration / (double) count;
        System.out.printf("Rate %g %n", rate);
    }
}
