/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import net.maizegenetics.gbs.maps.PETagsOnPhysicalMap;
import net.maizegenetics.gbs.maps.TagsOnGeneticMap;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMapV3;
import net.maizegenetics.gbs.tagdist.PETagCounts;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.util.MultiMemberGZIPInputStream;

/**
 *
 * @author Fei Lu
 */
public class FeiUtils {
    
    public FeiUtils () {
        
    }
    
    public void mergeBlastResult (String inputDirS, String outputFileS) {
        File[] files = new File(inputDirS).listFiles();
        try {
            
        }
        catch (Exception e) {
            
        }
    }
    
    public void splitFastaFileS (String inputFileS, String outputDir, int fileNum) {
        try{
            BufferedReader br;
            if (inputFileS.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(inputFileS)))));
            } else {
                br = new BufferedReader(new FileReader(new File(inputFileS)), 65536);
            }
            String temp;
            int tagNum = 0;
            while ((temp = br.readLine()) != null) {
                tagNum++;
            }
            br.close();
            tagNum = tagNum/2;
            int[] size = new int[fileNum];
            int left = tagNum % fileNum;
            if (left == 0) {
                for (int i = 0; i < size.length; i++) {
                    size[i] = tagNum / fileNum;
                }
            }
            else {
                int base = tagNum/fileNum;
                for (int i = 0; i < left; i++) {
                    size[i] = base+1;
                }
                for (int i = left; i < fileNum; i++) {
                    size[i] = base;
                }
            }
            if (inputFileS.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(inputFileS)))));
            } else {
                br = new BufferedReader(new FileReader(new File(inputFileS)), 65536);
            }
            for (int i = 0; i < fileNum; i++) {
                String fileName = outputDir+"tag"+this.get5DigitNumber(i)+".fasta";
                BufferedWriter bw = new BufferedWriter (new FileWriter(fileName), 65536);
                for (int j = 0; j < size[i]; j++) {
                    bw.write(br.readLine());
                    bw.newLine();
                    bw.write(br.readLine());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    private String get5DigitNumber (int j) {
        String numS = String.valueOf(j);
        int repeat = 5 - numS.length();
        for (int i = 0; i < repeat; i++) {
            numS = "0"+numS;
        }
        return numS;
    }
    
    public void convertSAM2TOPM (String infileS, String outfileS) {
        SAMConverterPlugin scp = new SAMConverterPlugin ();
        String arguments = "-i " + infileS + " -t -o " + outfileS;
        String[] args = arguments.split(" ");
		scp.setParameters(args);
		scp.performFunction(null);
    }
    
    /**
     * Convert TagCounts to Fastq for alignment
     * @param inputFileS
     * @param outputFileS 
     */
    public void convertTagCount2Fastq (String inputFileS, String outputFileS) {
        TagCountToFastqPlugin umithm = new TagCountToFastqPlugin();
		String arguments = "-i " + inputFileS + " -o " + outputFileS;
		String[] args = arguments.split(" ");
		umithm.setParameters(args);
		umithm.performFunction(null);       
    }
    
    public void convertTOPM2Fastq (String topmFileS, String outputFileS) {
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmFileS);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outputFileS), 65536);
            for (int i = 0; i < topm.getTagCount(); i++) {
                bw.write("@length="+topm.getTagLength(i)+"count="+1);
                bw.newLine();
                long[] tag = topm.getTag(i);
                bw.write(BaseEncoder.getSequenceFromLong(tag).substring(0, topm.getTagLength(i)));
                bw.newLine();
                bw.write("+");
                bw.newLine();
                for (int j = 0; j < topm.getTagLength(i); j++) {
                    bw.write("f");
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void mkSmallPETagCountsFile (String inputFileS, String outputFileS, int startIndex, int size) {
        PETagCounts ptc = new PETagCounts (inputFileS, FilePacking.Bit);
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileS), 65536));
            dos.writeInt(ptc.getTagSizeInLong());
            dos.writeInt(size);
            for (int i = startIndex; i < startIndex+size; i++) {
                
                for (int j = 0; j < ptc.getTagSizeInLong(); j++) {
                    dos.writeLong(ptc.getTagF(i)[j]);
                }
                dos.writeShort(ptc.getTagFLength(i));
                for (int j = 0; j < ptc.getTagSizeInLong(); j++) {
                    dos.writeLong(ptc.getTagB(i)[j]);
                }
                dos.writeShort(ptc.getTagBLength(i));
                dos.writeByte(ptc.getContigLengthInLong(i));
                for (int j = 0; j < ptc.getContigLengthInLong(i); j++) {
                    dos.writeLong(ptc.getContig(i)[j]);
                }
                dos.writeShort(ptc.getContigLength(i));
                dos.writeInt(ptc.getReadCount(i));
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Only output N row tagCounts
     * @param inputFileS
     * @param outputFileS
     * @param tagNum 
     */        
    public void mkSmallTagCountsFile (String inputFileS, String outputFileS, int startIndex, int tagNum) {
        TagCounts tc = new TagCounts (inputFileS, TagsByTaxa.FilePacking.Bit);
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileS), 65536));
            fw.writeInt(tagNum);
            fw.writeInt(tc.getTagSizeInLong());
            for (int i = startIndex; i < startIndex+tagNum; i++) {
                for (int j = 0; j < tc.getTagSizeInLong(); j++) {
                    fw.writeLong(tc.getTag(i)[j]);
                }
                fw.writeByte(tc.getTagLength(i));
                fw.writeInt(tc.getReadCount(i));
            }
            fw.close();
        } catch (Exception e) {
            System.out.println("Catch in writeTextDistFile writing output file e=" + e);
            e.printStackTrace();
        }
    }
    
    public void mkAlignmentCompareTable (String TOGMFileS, String topmFileS, String ptopmFileS, String outfileS) {
        int tagLimit = 30000;
        TagsOnGeneticMap togm = new TagsOnGeneticMap (TOGMFileS, TagsByTaxa.FilePacking.Text);
        TagsOnPhysicalMap topm = new TagsOnPhysicalMap (topmFileS, true); 
        PETagsOnPhysicalMap ptopm = new PETagsOnPhysicalMap (ptopmFileS, TagsByTaxa.FilePacking.Text);
        int cnt = 0;
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter (outfileS), 65536);
            bw.write("Tag\tGChr\tGPOs\tPEChr\tPEPos\tTChr\tTPos");
            bw.newLine();
            for (int i = 0; i < ptopm.getTagCount(); i++) {
                //if (ptopm.getChrF(i) == Integer.MIN_VALUE) continue;
                long[] PEtag = ptopm.getTagF(i);
                long[] tag = new long[togm.getTagSizeInLong()];
                for (int j = 0; j < tag.length; j++) {
                    tag[j] = PEtag[j];
                }
                int index1 = togm.getTagIndex(tag);
                if (index1 < 0) continue;
                int index2 = topm.getTagIndex(tag);
                if (index2 < 0) continue;
                bw.write(BaseEncoder.getSequenceFromLong(togm.getTag(i))+"\t");
                bw.write(String.valueOf(togm.getGChr(index1))+"\t"+String.valueOf(togm.getGPos(index1))+"\t");
                bw.write(String.valueOf(ptopm.getChrF(i))+"\t"+String.valueOf(ptopm.getPosStartF(i))+"\t");
                bw.write(String.valueOf(topm.getChromosome(index2))+"\t"+String.valueOf(topm.getPositionArray(index2)[2]));
                bw.newLine();
                cnt++;
                if (cnt> tagLimit) break;
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
}
