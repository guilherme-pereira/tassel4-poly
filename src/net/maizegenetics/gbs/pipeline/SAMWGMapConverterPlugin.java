/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.pipeline;
import java.awt.Frame;
import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.pal.alignment.BitAlignmentHDF5;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;
import org.apache.log4j.Logger;

/**
 * This class reads in SAM mapping results tests them against an anchor map 
 * and creates a update HDF5 TOPM file
 *
 * @author buckler
 *
 */
public final class SAMWGMapConverterPlugin extends AbstractPlugin{
    boolean cleanCutSites=true;
    private static final Logger myLogger=Logger.getLogger(SAMWGMapConverterPlugin.class);
    private static ArgsEngine myArgsEngine;
    private static String inputFileName=null;
    private static String outputFileName=null;
    private boolean textFormat = false;
    private int tagLengthInLong = 2;
    private BitAlignmentHDF5[] anchAlign;

    public SAMWGMapConverterPlugin(){
        super(null, false);
    }

    public SAMWGMapConverterPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        myLogger.info(
           "\n\nUsage is as follows:\n"
            + "-i  Name of input file in SAM text format (required)\n"
            + "-a  Name of anchor maps in HDF5 format"
            + "-o  Name of output file (default output.topm.bin)\n"
            + "-t  Specifies text output format\n"
            + "-l  tag length in integer multiples of 32 bases (default=2)\n\n"    
        );
    }

    @Override
    public DataSet performFunction(DataSet input) {
        TagsOnPhysicalMap topm = new TagsOnPhysicalMap();
        topm.readSAMFile(inputFileName, tagLengthInLong);
        topm.sort();
        try {
            if(textFormat==true){
                topm.writeTextFile(new File(outputFileName));
            } else {
                topm.writeBinaryFile(new File(outputFileName));
            }
        } catch (Exception e) {
            System.out.println("Catch in writing binary topm file: " + e);
        }
        writeLogFile(topm);
        return null;
    }
    
    private void writeLogFile(TagsOnPhysicalMap topm) {
        try {
            DataOutputStream report = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileName+".log"), 65536));
            int[] aligned=topm.mappedTags();
            int unique=0, multi=1;  // the indices of aligned
            int unaligned=topm.getTagCount()-aligned[unique]-aligned[multi];
            report.writeBytes(
                "Input file: "+inputFileName+"\n"+
                "Output file: "+outputFileName+"\n"+
                "Total "+topm.getTagCount()+" tags\n\t"
                    +aligned[unique]+" were aligned to unique postions\n\t"
                    +aligned[multi]+" were aligned to multiple postions\n\t"
                    +unaligned+" could not be aligned.\n\n"
            );
            int[] dist = topm.mappingDistribution();
            report.writeBytes("nPositions  nTags\n");
            for (int i = 0; i < dist.length; i++) {
                if (dist[i]>0) {
                    if (i<10) {
                        report.writeBytes(i+"           "+dist[i]+"\n");
                    } else if (i<100) {
                        report.writeBytes(i+"          "+dist[i]+"\n");
                    } else if (i<1000) {
                        report.writeBytes(i+"         "+dist[i]+"\n");
                    }
                }
            }
            report.close();
        } catch(Exception e) {
            myLogger.warn("Caught exception while writing log file: "+e);
        }
    }

    @Override
    public void setParameters(String[] args) {
        if(args.length==0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-i", "--input-file", true);
            myArgsEngine.add("-o", "--output-file", true);
            myArgsEngine.add("-t", "--text-format");
            myArgsEngine.add("-l", "--tag-length-in-mutiples-of-32-bases", true);
        }
        myArgsEngine.parse(args);

        if(myArgsEngine.getBoolean("-i")){
            File inputFile = new File(myArgsEngine.getString("-i"));
            if (!inputFile.isFile()) {
                printUsage();
                throw new IllegalArgumentException("The input name you supplied is not a valid file.");
            }
            inputFileName = inputFile.getAbsolutePath();
            outputFileName = inputFile.getParent()+File.pathSeparator+"output.topm.bin";
        }else{
            printUsage();
            throw new IllegalArgumentException("Please supply an input file name.");
        }
        if(myArgsEngine.getBoolean("-t")) {
            textFormat=true;
            outputFileName = outputFileName.replace(".bin", ".txt");
        }
        if(myArgsEngine.getBoolean("-o")){
            outputFileName = myArgsEngine.getString("-o");
        }
        if(myArgsEngine.getBoolean("-l")){
            tagLengthInLong = Integer.parseInt(myArgsEngine.getString("-l"));
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
    
    public static void main(String[] args) {
        String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/";
        int sC=8;
        int eC=10;
        String hmpfile=root+"temp_Ed/July_2012_Build.BPEC.Highf.c#.hmp.txt";
        String h5file=root+"test/July_2012_Build.BPEC.Highf.c#.hmp.h5";
//        for (int i = sC; i <= eC; i++) {
//            String infile=hmpfile.replace("#", ""+i);
//            SBitAlignment a=(SBitAlignment)ImportUtils.readFromHapmap(infile, null);
//            String h5out=h5file.replace("#", ""+i);
//            SBitAlignmentNucleotideHDF5.createFile(a, h5out);
//        }
//        System.exit(0);
        
//        SBitAlignmentNucleotideHDF5[] anchAlign=new SBitAlignmentNucleotideHDF5[13];
//        for (int i = sC; i <= eC; i++) {
//            String h5out=h5file.replace("#", ""+i);
//            anchAlign[i]=SBitAlignmentNucleotideHDF5.getInstance(h5out,""+i);
//            System.out.println("chr"+i+": sites="+anchAlign[i].getSiteCount());
//        }
//        System.exit(0);
       
        String samFile=root+"04_TOPM/Bowtie2/h100kAllZeaMasterTags_c10_20120626_k.sam";
        String topm=root+"test/topm100.txt";
        
        args = new String[] {
            "-i", samFile,
            "-o", topm,
       //     "-a", h5file,
            "-t",
//            "-sC", ""+sC,  // Start chromosome
//            "-eC", ""+eC // End chromosome
        };
        SAMWGMapConverterPlugin plugin = new SAMWGMapConverterPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }


}
