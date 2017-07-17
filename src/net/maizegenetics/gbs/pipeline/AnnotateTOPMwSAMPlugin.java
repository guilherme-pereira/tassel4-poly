/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

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
import javax.swing.ImageIcon;
import net.maizegenetics.annotation.UnderConstruction;
import net.maizegenetics.gbs.maps.TagMappingInfo;
import net.maizegenetics.gbs.maps.TagsOnPhysMapHDF5;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.MultiMemberGZIPInputStream;
import org.apache.log4j.Logger;

/**
 * This class reads in SAM mapping results tests them against an anchor map 
 * and creates a updated HDF5 TOPM file.
 *
 * TODO:
 * <li> Add mapping information from Bowtie2
 * <li> Add mapping information from BWA
 * <li> Add mapping information from BLAST?
 * <li> Run genetic to compare hypotheses
 * <li> Call resort
 * @author Ed Buckler and Fei Lu
 *
 */
@UnderConstruction(owner="Fei Lu", otherContacts="Ed Buckler")
public final class AnnotateTOPMwSAMPlugin extends AbstractPlugin{
    boolean cleanCutSites=true;
    private static final Logger myLogger=Logger.getLogger(AnnotateTOPMwSAMPlugin.class);
    private static ArgsEngine myArgsEngine;
    private static String inputFileName=null;
    private static String topmFileName=null;
    private boolean textFormat = false;
    private int tagLengthInLong = 2;
   // private SBitAlignmentNucleotideHDF5[] anchAlign;

    public AnnotateTOPMwSAMPlugin(){
        super(null, false);
        
    }

    public AnnotateTOPMwSAMPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        myLogger.info(
           "\n\nUsage is as follows:\n"
            + "-i  Name of input file in SAM text format (required)\n"
            + "-a  Name of anchor maps in HDF5 format"
            + "-m  Name of topm hdf5 file (default output.topm.bin)\n"
            + "-t  Specifies text output format\n"
            + "-l  tag length in integer multiples of 32 bases (default=2)\n\n"    
        );
    }

    @Override
    public DataSet performFunction(DataSet input) {
        TagsOnPhysMapHDF5 topm = new TagsOnPhysMapHDF5(topmFileName,true);
        for (int i = 0; i < topm.getTagCount(); i++) {
            topm.setMultimaps(i, (byte)0);
         //   System.out.printf("%d %d %n",i,topm.getMultiMaps(i));
            
        }
        readSAMFile(topm, inputFileName, tagLengthInLong);
        topm.sort();
        
        //writeLogFile(topm);
        return null;
    }
    
    /** Reads SAM files output from bowtie2 */
    public void readSAMFile(TagsOnPhysMapHDF5 topm, String inputFileName, int tagLengthInLong) {
        System.out.println("Reading SAM format tag alignment from: " + inputFileName);
        this.tagLengthInLong = tagLengthInLong;
        String inputStr = "Nothing has been read from the file yet";
       // int nHeaderLines = countTagsInSAMfile(inputFileName); // detects if the file is bowtie2, initializes topm matrices
        int tagIndex = Integer.MIN_VALUE;
        try {
            BufferedReader br;
            if (inputFileName.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(inputFileName)))));
            } else {
                br = new BufferedReader(new FileReader(new File(inputFileName)), 65536);
            }
            while(br.readLine().startsWith("@PG")==false) {};
            int count=0;
            while((inputStr = br.readLine())!=null) {
                String[] d=inputStr.split("\\s");
                int orientiation=Integer.parseInt(d[1]);
                if(orientiation==4) continue; //no alignment
                int chr=Integer.parseInt(d[2]);
                int position=Integer.parseInt(d[3]);
                byte strand=1;
                String seqS=d[9];
                if((orientiation==16)||(orientiation==272)) {
                    seqS=BaseEncoder.getReverseComplement(seqS);
                    strand=-1;
                    //may need to change position...
                }
                long[] seq=BaseEncoder.getLongArrayFromSeq(seqS,tagLengthInLong*32);
                int index=topm.getTagIndex(seq);
                if(index<0) {System.out.println(count+"Tag not found:"+inputStr);}
     //           else {System.out.println(count+"Tag is  found:"+inputStr);}
                byte alignScore=(byte)(128-Integer.parseInt(d[11].split(":")[2]));
                count++;
                System.out.println(count+"Tag is  found:"+inputStr);
                TagMappingInfo theTMI=new TagMappingInfo(chr,strand,position, position+64,alignScore);
                int mm=topm.getMultiMaps(index)+1;
                System.out.printf("TagIndex %d  MultiMap: %d %n",index,mm);
                topm.setAlternateTagMappingInfo(index, mm-1, theTMI);
                
             //   System.out.println(count++);
                
            }
            br.close();
            topm.getFileReadyForClosing();
        } catch (Exception e) {
            System.out.println("\n\nCatch in reading SAM alignment file at tag " + tagIndex + ":\n\t" + inputStr + "\nError: " + e + "\n\n");
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void writeLogFile(TagsOnPhysicalMap topm) {
        try {
            DataOutputStream report = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(topmFileName+".log"), 65536));
            int[] aligned=topm.mappedTags();
            int unique=0, multi=1;  // the indices of aligned
            int unaligned=topm.getTagCount()-aligned[unique]-aligned[multi];
            report.writeBytes(
                "Input file: "+inputFileName+"\n"+
                "Output file: "+topmFileName+"\n"+
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
            myArgsEngine.add("-m", "--topm-file", true);
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
            topmFileName = inputFile.getParent()+File.pathSeparator+"output.topm.bin";
        }else{
            printUsage();
            throw new IllegalArgumentException("Please supply an input file name.");
        }
        if(myArgsEngine.getBoolean("-t")) {
            textFormat=true;
            topmFileName = topmFileName.replace(".bin", ".txt");
        }
        if(myArgsEngine.getBoolean("-m")){
            topmFileName = myArgsEngine.getString("-m");
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
//        String inTOPMFile=root+"04_TOPM/Bowtie2/AllZeaMasterTags_c10_20120703.topm";
//        boolean binary=!inTOPMFile.endsWith(".txt");
//        TagsOnPhysicalMap inTOPM = new TagsOnPhysicalMap(inTOPMFile, binary);
//        TagsOnPhysMapHDF5.createFile(inTOPM, root+"04_TOPM/Bowtie2/AllZeaMasterTags_c10_20120703.topm.h5",4,8);
//        inTOPM.printRows(100000);

        
        String topm=root+"04_TOPM/Bowtie2/AllZeaMasterTags_c10_20120703.topm.h5";
        
        args = new String[] {
            "-i", samFile,
            "-m", topm,
       //     "-a", h5file,
            "-t",
//            "-sC", ""+sC,  // Start chromosome
//            "-eC", ""+eC // End chromosome
        };
        AnnotateTOPMwSAMPlugin plugin = new AnnotateTOPMwSAMPlugin();
        plugin.setParameters(args);
        plugin.performFunction(null);
    }


}
