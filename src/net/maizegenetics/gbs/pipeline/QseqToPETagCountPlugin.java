package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.BufferedOutputStream;

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Arrays;

import net.maizegenetics.util.MultiMemberGZIPInputStream;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.homology.PEParseBarcodeRead;
import net.maizegenetics.gbs.homology.PEReadBarcodeResult;
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
 * Derives a PETagCount list for a pair of Qseq files. The forward and backward tags are ordered during processing
 *
 * Keeps only good reads having a barcode and a cut site and no N's in the useful part of the sequence. 
 * For the barcoded end, trims off the barcodes and truncates sequences
 * that (1) have a second cut site, or (2) read into the common adapter.
 * For the unbarcoded end, trims off the barcodes and truncates sequences
 * that (1) have a second cut site, or (2) read into the barcode adapter.
 */
public class QseqToPETagCountPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(QseqToPETagCountPlugin.class);    
    String inputFileSWithBarcode = null;
    String inputFileSWithoutBarcode = null;
    String keyfile = null;
    String enzyme = null;
    int tagLengthInLong = 8;
    int minCount = 1;
    String outputDirS = null;

    public QseqToPETagCountPlugin() {
        super(null, false);
    }

    public QseqToPETagCountPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -iF Qseq file with barcode\n"
                + " -iB Qseq file without barcode\n"
                + " -k  Key file listing barcodes for each sample\n"
                + " -e  Enzyme used to create the GBS library, if it differs from the one listed in the key file.\n"
                + " -l  Tag length in Long type (Default is 2).\n"
                + " -c  Minimum tag count (default is 1).\n"
                + " -o  Output directory to contain .pe.cnt files\n");
    }

    @Override
    public DataSet performFunction(DataSet input) {
        countTags(keyfile, enzyme, inputFileSWithBarcode, inputFileSWithoutBarcode, minCount, tagLengthInLong);  //TODO change to perform function
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
            engine.add("-iF", "--input-barcode", true);
            engine.add("-iB", "--input-noBarcode", true);
            engine.add("-k", "--key-file", true);
            engine.add("-e", "--enzyme", true);
            engine.add("-l", "--tag-lengthInLong", true);
            engine.add("-c", "--min-count", true);
            engine.add("-o", "--output-directory", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-iF")) {
            inputFileSWithBarcode = engine.getString("-iF");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the Illumina file with barcode");
        }

        if (engine.getBoolean("-iB")) {
            inputFileSWithoutBarcode = engine.getString("-iB");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the Illumina file without barcode");
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
        }

        if (engine.getBoolean("-l")) {
            tagLengthInLong = Integer.parseInt(engine.getString("-l"));
        }
         
        if (engine.getBoolean("-c")) {
            minCount = Integer.parseInt(engine.getString("-c"));
        }

        if (engine.getBoolean("-o")) {
            outputDirS = engine.getString("-o");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a PETagCount directory.");
        }
    }

    /**
     * Derives a PETagCount list for a pair of Qseq files.
     *
     * @param keyFileS  A key file (a sample key by barcode, with a plate map included).
     * @param enzyme  The enzyme used to create the library (currently ApeKI or PstI).
     * @param inputFileSWithBarcode  Filename of input Qseq file with barcode
     * @param inputFileSWithoutBarcode  Filename of input Qseq file without barcode
     * @param minCount  The minimum number of occurrences of a tag in a Illumina file for it to be included in the output tagCounts file
     * @param tagLengthInLong  Tag Length in Long primitive data type, the default is 8 (8*32=256 bp)
     */
    public void countTags(String keyFileS, String enzyme, String inputFileSWithBarcode, String inputFileSWithoutBarcode, int minCount, int tagLengthInLong) {
        BufferedReader brF;
        BufferedReader brB;
        String[] countFileNames = null;
        int allReads = 0, goodBarcodedReads = 0;
        File inputFileWithBarcode = new File(inputFileSWithBarcode);
        File inputFileWithoutBarcode = new File(inputFileSWithoutBarcode);
        String[] filenameField = inputFileWithBarcode.getName().split("_");
        PEParseBarcodeRead thePBR;
        if(filenameField.length==4) {thePBR=new PEParseBarcodeRead(keyFileS, enzyme, filenameField[0], filenameField[1]);}
        else {
            System.out.println("Error in parsing file name:");
            System.out.println("The filename does not contain either 3 or 5 underscore-delimited values.");
            System.out.println("Expect: flowcell_lane_qseq.txt OR code_flowcell_s_lane_qseq.txt");
            System.out.println("Filename: "+inputFileSWithBarcode);
            return;
        }
        System.out.println("Total barcodes found in lane:"+thePBR.getBarCodeCount());
        if(thePBR.getBarCodeCount() == 0){
            System.out.println("No barcodes found.  Skipping this flowcell lane.");
            System.exit(1);
        }
        String[] taxaNames=new String[thePBR.getBarCodeCount()];
        for (int i = 0; i < taxaNames.length; i++) {
            taxaNames[i]=thePBR.getTheBarcodes(i).getTaxaName();
        }
		Arrays.sort(taxaNames);
        
        TaxonOutput[] tops = new TaxonOutput[taxaNames.length];
        File outputDir = new File(outputDirS);
        for (int i = 0; i < tops.length; i++) {
            tops[i] = new TaxonOutput(outputDir, taxaNames[i], tagLengthInLong);
        }
        
        try {
            if (inputFileWithBarcode.getName().endsWith(".gz")) {
                brF = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(inputFileWithBarcode))));
                brB = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(inputFileWithoutBarcode))));
            } else {
                brF = new BufferedReader(new FileReader(inputFileWithBarcode), 65536);
                brB = new BufferedReader(new FileReader(inputFileWithoutBarcode), 65536);
            }
            String sequenceF = null, qualityScoreF = null, temp = null;
            String sequenceB = null, qualityScoreB = null;
            while (((temp = brF.readLine()) != null)) {
                String[] temF = temp.split("\t");
                String[] temB = brB.readLine().split("\t");
                allReads++;
                try {
                    sequenceF = temF[8];
                    qualityScoreF = temF[9];
                    sequenceB = temB[8];
                    qualityScoreB = temB[9];
                } catch (NullPointerException e) {
                    System.out.println("Read a line that lacks a sequence and "
                            + "quality score in fields 9 and 10.  Your file may have been corrupted.");
                    System.exit(0);
                }
                PEReadBarcodeResult rr = thePBR.parseReadIntoTagAndTaxa(sequenceF, qualityScoreF, sequenceB,qualityScoreB, false, 0, tagLengthInLong);
                if (rr != null) {
                    goodBarcodedReads++;
                    int index = Arrays.binarySearch(taxaNames, rr.getTaxonName());
                    tops[index].write(rr);
                }
                if (allReads % 1000000 == 0) {
                    System.out.println("Total Reads:" + allReads + " Reads with barcode and cut site overhang:" + goodBarcodedReads);
                 }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        for (int i = 0; i < tops.length; i++) {
            tops[i].close();
        }
    }
    
    private class TaxonOutput {
        File taxonoutfile;
        DataOutputStream dos;
        
        TaxonOutput (File parentDir, String taxonName, int tagLengthInLong) {
            setupOutput(parentDir, taxonName, tagLengthInLong);
        }
        
        void setupOutput (File parentDir, String taxonName, int tagLengthInLong) {
            taxonoutfile = new File (parentDir, taxonName.replaceAll(":", "_")+".pe.cnt");
            try {
                dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(taxonoutfile), 65536));
                dos.writeInt(tagLengthInLong);
                dos.writeInt(-1); // -1 means parsed read without sorting, merging and contiging, tagNum is not known neither
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        
        void write (PEReadBarcodeResult rr) {
            try {
                for (int i = 0; i < rr.getTagLengthInLong(); i++) {
                    dos.writeLong(rr.getReadF()[i]);
                }
                dos.writeShort(rr.getLengthF());
                for (int i = 0; i < rr.getTagLengthInLong(); i++) {
                    dos.writeLong(rr.getReadB()[i]);
                }
                dos.writeShort(rr.getLengthB());
                dos.writeByte(0); // 0 means contig doesn't exist.
                dos.writeShort(0);
                dos.writeInt(1);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        
        void close () {
            try {
                dos.flush();
                dos.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
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
