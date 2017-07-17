/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.maps;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import net.maizegenetics.gbs.tagdist.AbstractPETags;
import net.maizegenetics.gbs.tagdist.PETagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.gbs.util.SAMUtils;

/**
 *
 * @author Fei Lu
 * @deprecated 
 */
public class PETagsOnPhysicalMap extends AbstractPETagsOnPhysicalMap {
    
    
    public PETagsOnPhysicalMap (String infileS, FilePacking format) {
        this.readDistFile(infileS, format);
    }
    
    public PETagsOnPhysicalMap (PETagCounts ptc, String fSamFileS, String bSamFileS, String contigSamFileS) {
        this.iniMatrix(ptc.getTagSizeInLong(), ptc.getTagCount());
        this.initializeValue(ptc);
        this.readSamFile(fSamFileS, PETagType.Forward);
        this.readSamFile(bSamFileS, PETagType.Backward);
        this.readSamFile(contigSamFileS, PETagType.Contig);
    }
    
    private void readSamFile (String samFileS, PETagType t) {
        int[] chromosome;
        int[] posStart;
        int[] posEnd;
        byte[] strand;
        if (t == PETagType.Forward) {
            chromosome = this.chrF;
            posStart = this.posStartF;
            posEnd = this.posEndF;
            strand = this.strandF;
        }
        else if (t == PETagType.Backward) {
            chromosome = this.chrB;
            posStart = this.posStartB;
            posEnd = this.posEndB;
            strand = this.strandB;
        }
        else {
            chromosome = this.chrContig;
            posStart = this.posStartContig;
            posEnd = this.posEndContig;
            strand = this.strandContig;
        }
        try {
            BufferedReader br = new BufferedReader(new FileReader(samFileS), 65536);
            String temp;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("@")) continue;
                this.parseSAMAlignment(temp, chromosome, posStart, posEnd, strand);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void parseSAMAlignment (String inputStr, int[] chromosome, int[] posStart, int[] posEnd, byte[] strand) {
        String[] temp = inputStr.split("\t");
        int tagIndex = Integer.parseInt(temp[0]);
        if (temp[2].equals("*")) return;
        HashMap<String, Integer> SAMFields = parseOptionalFieldsFromSAMAlignment(temp);
        byte bestHits = (byte) Math.min(SAMFields.get("nBestHits"), Byte.MAX_VALUE);
        if (bestHits != 1) return;
        chromosome[tagIndex] = Integer.parseInt(temp[2]);
        strand[tagIndex] = (Integer.parseInt(temp[1]) == 16) ? (byte) -1 : (byte) 1;
        this.recordStartEndPostionFromSAMAlign(tagIndex, strand[tagIndex], Integer.parseInt(temp[3]), temp[5], posStart, posEnd);
    }
    
    private void recordStartEndPostionFromSAMAlign(int tagIndex, byte strand, int pos, String cigar, int[] posStart, int[] posEnd) {
        int[] alignSpan = SAMUtils.adjustCoordinates(cigar, pos);
        try {
            if (strand == 1) {
                posStart[tagIndex] = alignSpan[0];
                posEnd[tagIndex] = alignSpan[1];
            } else if (strand == -1) {
                posStart[tagIndex] = alignSpan[1];
                posEnd[tagIndex] = alignSpan[0];
            } else {
                throw new Exception("Unexpected value for strand: " + strand + "(expect 1 or -1)");
            }
        } catch (Exception e) {
            System.out.println("Error in recordStartEndPostionFromSAMAlign: " + e);
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private HashMap<String, Integer> parseOptionalFieldsFromSAMAlignment(String[] inputLine) {
        HashMap<String, Integer> SAMFields = new HashMap<String, Integer>();
        for (int field = 11; field < inputLine.length; field++) { // Loop through all the optional field of the SAM alignment
            if (inputLine[field].regionMatches(0, "AS", 0, 2)) {        // AS = SAM format for alignment score of the best alignment.  Specific to bowtie2.
                SAMFields.put("bestScore", Integer.parseInt(inputLine[field].split(":")[2]));
            } else if (inputLine[field].regionMatches(0, "XS", 0, 2)) { // XS = SAM format for alignment score of 2nd best alignment.  Specific to bowtie2.
                SAMFields.put("nextScore", Integer.parseInt(inputLine[field].split(":")[2]));
            } else if (inputLine[field].regionMatches(0, "NM", 0, 2)) { // NM = SAM format for edit distance to the reference.  Common to BWA and Bowtie2.
                SAMFields.put("editDist", Integer.parseInt(inputLine[field].split(":")[2]));
            }
        }
        if (SAMFields.containsKey("bestScore")) {
            if (SAMFields.containsKey("nextScore")) {
                if (SAMFields.get("bestScore") > SAMFields.get("nextScore")) {
                     SAMFields.put("nBestHits", 1);
                } else {
                    SAMFields.put("nBestHits", (int)Byte.MAX_VALUE);  // Byte.MAX_VALUE will stand for an unknown # of multiple hits
                }
            } else {
                SAMFields.put("nBestHits", 1);
            }
        }
        return SAMFields;
    }
    

    private void initializeValue (PETagCounts ptc) {
        tagLengthInLong = ptc.getTagSizeInLong(); 
        for (int i = 0; i < this.getTagCount(); i++) {
            tagsF[i] = ptc.getTagF(i);
            tagsB[i] = ptc.getTagB(i);
            tagFLength[i] = ptc.getTagFLength(i);
            tagBLength[i] = ptc.getTagBLength(i);
            contigLengthInLong[i] = ptc.getContigLengthInLong(i);
            contig[i] = ptc.getContig(i);
            contigLength[i] = ptc.getContigLength(i);
            chrF[i]= Integer.MIN_VALUE;
            chrB[i]= Integer.MIN_VALUE;
            chrContig[i]= Integer.MIN_VALUE;
            posStartF[i] = Integer.MIN_VALUE;
            posStartB[i] = Integer.MIN_VALUE;
            posStartContig[i] = Integer.MIN_VALUE;
            posEndF[i] = Integer.MIN_VALUE;
            posEndB[i] = Integer.MIN_VALUE;
            posEndContig[i] = Integer.MIN_VALUE;
            strandF[i] = Byte.MIN_VALUE;
            strandB[i] = Byte.MIN_VALUE;
            strandContig[i] = Byte.MIN_VALUE;
        }
    }
    
    public void readDistFile (String infileS, FilePacking format) {
        System.out.println("Reading PTOPM file to " + infileS);
        File infile = new File (infileS);
        switch (format) {
            case Text:
                readTextPTOPMFile(infile);
                break;
            default:
                readBinaryPTOPMFile(infile);
                break;
        }
        System.out.println("PTOPM file read. Tatol: " + this.getTagCount() + " PETags");
    }
    
    private void readBinaryPTOPMFile (File infile) {
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infile), 65536));
            tagLengthInLong = dis.readInt();
            int tagNum = dis.readInt();
            this.iniMatrix(tagLengthInLong, tagNum);
            for (int i = 0; i < tagNum; i++) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    tagsF[i][j] = dis.readLong();
                }
                tagFLength[i] = dis.readShort();
                chrF[i] = dis.readInt();
                posStartF[i] = dis.readInt();
                posEndF[i] = dis.readInt();
                strandF[i] = dis.readByte();
                for (int j = 0; j < tagLengthInLong; j++) {
                    tagsB[i][j] = dis.readLong();
                }
                tagBLength[i] = dis.readShort();
                chrB[i] = dis.readInt();
                posStartB[i] = dis.readInt();
                posEndB[i] = dis.readInt();
                strandB[i] = dis.readByte();
                contigLengthInLong[i] = dis.readByte();
                contig[i] = new long[contigLengthInLong[i]];
                for (int j = 0; j < contig[i].length; j++) {
                    contig[i][j] = dis.readLong();
                }   
                contigLength[i] = dis.readShort();
                chrContig[i] = dis.readInt();
                posStartContig[i] = dis.readInt();
                posEndContig[i] = dis.readInt();
                strandContig[i] = dis.readByte();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void readTextPTOPMFile (File infile) {
         try {
            BufferedReader br = new BufferedReader (new FileReader(infile), 65536);
            tagLengthInLong = Integer.valueOf(br.readLine());
            int tagNum = Integer.valueOf(br.readLine()); 
            this.iniMatrix(tagLengthInLong, tagNum);
            br.readLine();
            for (int i = 0; i < tagNum; i++) {
                String[] temp = br.readLine().split("\t");
                long[] t = BaseEncoder.getLongArrayFromSeq(temp[0]);
                for (int j = 0; j < tagLengthInLong; j++) {
                    tagsF[i][j] = t[j];
                }
                tagFLength[i] = Short.valueOf(temp[1]);
                chrF[i] = Integer.valueOf(temp[2]);
                posStartF[i] = Integer.valueOf(temp[3]);
                posEndF[i] = Integer.valueOf(temp[4]);
                strandF[i] = Byte.valueOf(temp[5]);
                t = BaseEncoder.getLongArrayFromSeq(temp[6]);
                for (int j = 0; j < tagLengthInLong; j++) {
                    tagsB[i][j] = t[j];
                }
                tagBLength[i] = Short.valueOf(temp[7]);
                chrB[i] = Integer.valueOf(temp[8]);
                posStartB[i] = Integer.valueOf(temp[9]);
                posEndB[i] = Integer.valueOf(temp[10]);
                strandB[i] = Byte.valueOf(temp[11]);
                contigLengthInLong[i] = Byte.valueOf(temp[12]);
                contigLength[i] = Short.valueOf(temp[13]);
                chrContig[i] = Integer.valueOf(temp[14]);
                posStartContig[i] = Integer.valueOf(temp[15]);
                posEndContig[i] = Integer.valueOf(temp[16]);
                strandContig[i] = Byte.valueOf(temp[17]);
                this.contig[i] = new long[contigLengthInLong[i]];
                if (contigLengthInLong[i] != 0) {
                    t = BaseEncoder.getLongArrayFromSeq(temp[18]);
                    for (int j = 0; j < t.length; j++) {
                        contig[i][j] = t[j];
                    }
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
   
    public void writeDistFile (String outfileS, FilePacking format) {
        System.out.println("Writing PTOPM file to " + outfileS);
        switch (format) {
            case Text:
                writeTextPTOPMFile(outfileS);
                break;
            default:
                writeBinaryPTOPMFile(outfileS);
                break;
        }
        System.out.println("PPTOM file written");
    }
    
    private void writeBinaryPTOPMFile (String outfileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
            dos.writeInt(tagLengthInLong);
            dos.writeInt(this.getTagCount());
            for (int i = 0; i < this.getTagCount(); i++) {
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    dos.writeLong(this.tagsF[i][j]);
                }
                dos.writeShort(this.tagFLength[i]);
                dos.writeInt(this.chrF[i]);
                dos.writeInt(this.posStartF[i]);
                dos.writeInt(this.posEndF[i]);
                dos.writeByte(this.strandF[i]);
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    dos.writeLong(this.tagsB[i][j]);
                }
                dos.writeShort(this.tagBLength[i]);
                dos.writeInt(this.chrB[i]);
                dos.writeInt(this.posStartB[i]);
                dos.writeInt(this.posEndB[i]);
                dos.writeByte(this.strandB[i]);
                dos.writeByte(this.contigLengthInLong[i]);
                for (int j = 0; j < this.contigLengthInLong[i]; j++) {
                    dos.writeLong(contig[i][j]);
                }
                dos.writeShort(this.contigLength[i]);
                dos.writeInt(this.chrContig[i]);
                dos.writeInt(this.posStartContig[i]);
                dos.writeInt(this.posEndContig[i]);
                dos.writeByte(this.strandContig[i]);
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void writeTextPTOPMFile (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write(String.valueOf(tagLengthInLong));
            bw.newLine();
            bw.write(String.valueOf(this.getTagCount()));
            bw.newLine();
            bw.write("TagF\tLengthF\tChrF\tPosStartF\tPosEndF\tStrandF\tTagB\tLengthB\tChrB\tPosStartB\tPosEndB\tStrandB\tLengthInLongContig\tLengthContig\tChrContig\tPosStartContig\tPosEndContig\tStrandContig\tTagContig");
            bw.newLine();
            for (int i = 0; i < this.getTagCount(); i++) {
                bw.write(BaseEncoder.getSequenceFromLong(this.getTagF(i))+"\t"+String.valueOf(this.getTagFLength(i))+"\t");
                bw.write(String.valueOf(chrF[i])+"\t"+String.valueOf(posStartF[i])+"\t"+String.valueOf(posEndF[i])+"\t"+String.valueOf(strandF[i])+"\t");
                bw.write(BaseEncoder.getSequenceFromLong(this.getTagB(i))+"\t"+String.valueOf(this.getTagBLength(i))+"\t");
                bw.write(String.valueOf(chrB[i])+"\t"+String.valueOf(posStartB[i])+"\t"+String.valueOf(posEndB[i])+"\t"+String.valueOf(strandB[i])+"\t");
                bw.write(String.valueOf(this.contigLengthInLong[i])+"\t"+String.valueOf(this.contigLength[i])+"\t");
                bw.write(String.valueOf(chrContig[i])+"\t"+String.valueOf(posStartContig[i])+"\t"+String.valueOf(posEndContig[i])+"\t"+String.valueOf(strandContig[i])+"\t");
                if (this.getContigLengthInLong(i) != 0) bw.write(BaseEncoder.getSequenceFromLong(this.contig[i])+"\t");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
