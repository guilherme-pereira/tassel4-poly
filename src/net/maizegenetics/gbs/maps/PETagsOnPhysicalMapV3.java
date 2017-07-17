/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.maps;

import cern.colt.GenericSorting;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.gbs.util.SAMUtils;
import net.maizegenetics.util.MultiMemberGZIPInputStream;

/**
 * This class hold the multiple mapping positions of PE tags, including forward and backward PE tags. Data in this class are used to annotate HDF5TOPM
 * When contig exist, the forward and backword tags are the contig and reverse complement of the contig, respectively.
 * The PE tags were first aligned with bowtie2 (-k N), any number of mapping position can be imported
 * PE tags of both end are stored in one long array with variable length, then they were truncated to 2 longs. This is easy for searching
 * Using a pairIndex, each 64 bp tag has a index pointing to the other end of PE.
 * It is worth noting that, the 64 bp tags are not unique. Rarely, there are multiple identical 64 bp tags, although the full length PE tags are unique.
 * When multiple identical 64 bp tags exist, the one with longest full length PE is used to annotate HDF5TOPM.
 * @author Fei Lu
 */
public class PETagsOnPhysicalMapV3 implements Tags {
    protected int tagLengthInLong = 2;
    long[][] tags;
    byte[] variableTagLengthInLong;
    short[] tagLength;
    boolean[] ifContig;
    int[] pairIndex;
    byte[] mappingNum;
    byte[][] chr;
    byte[][] strand;
    int[][] startPos;
    short[][] score;
    byte[][] divergence;
    
    /**
     * Constructor from a file
     * @param PETOPMFileS 
     */
    public PETagsOnPhysicalMapV3 (String PETOPMFileS) {
        this.readBinaryFile(PETOPMFileS);
    }
    
    /**
     * Constructor using fasta files and sam files
     * @param fFastaFileS
     * @param bFastaFileS
     * @param fSamFileS
     * @param bSamFileS 
     */
    public PETagsOnPhysicalMapV3 (String fFastaFileS, String bFastaFileS, String fSamFileS, String bSamFileS) {
        int tagNum = 0;
        try {
            BufferedReader br = new BufferedReader (new FileReader(fFastaFileS), 65536);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                tagNum++;
            }
            System.out.println(String.valueOf(tagNum)+" PE tags (one end) in total");
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.iniMatrix(tagNum);
        this.importSamFile(fFastaFileS, fSamFileS, 0);
        this.importSamFile(bFastaFileS, bSamFileS, tagNum/2);
        this.trancateTag();
        sort();
/*        
        int tagIndex = 60;
        System.out.println(BaseEncoder.getSequenceFromLong(tags[tagIndex]));
        System.out.println(tagLength[tagIndex]+"\t"+chr[tagIndex][0]+"\t"+startPos[tagIndex][0]);
        System.out.println(BaseEncoder.getSequenceFromLong(tags[pairIndex[tagIndex]]));
        System.out.println(tagLength[pairIndex[tagIndex]]+"\t"+chr[pairIndex[tagIndex]][0]+"\t"+startPos[pairIndex[tagIndex]][0]);
        System.out.println("Good");
        String a = BaseEncoder.getSequenceFromLong(tags[tagIndex]);
        long[] t = BaseEncoder.getLongArrayFromSeq(a);
        System.out.println(this.getTagIndex(t));
*/
    }
    
    /**
     * Full length PE tags are truncated to certain size (2 longs in GBS Build V3.0)
     */
    public void trancateTag () {
        String nullS = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        for (int i = 0; i < tags.length; i++) {
            long[] t = new long[tagLengthInLong];
            if (variableTagLengthInLong[i] < tagLengthInLong) {
                for (int j = 0; j < variableTagLengthInLong[i]; j++) {
                    t[j] = tags[i][j];
                }
                for (int j = variableTagLengthInLong[i]; j < tagLengthInLong; j++) {
                    t[j] = BaseEncoder.getLongFromSeq(nullS);
                }
            }
            else {
                for (int j = 0; j < tagLengthInLong; j++) {
                    t[j] = tags[i][j];
                }
            }
            tags[i] = t;
        }
        System.out.println("PE tags were truncated to " + String.valueOf(tagLengthInLong*BaseEncoder.chunkSize) + " bp");
    }
    
    /**
     * Write to a file
     * @param outfileS 
     */
    public void writeBinaryFile (String outfileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
            dos.writeInt(tagLengthInLong);
            dos.writeInt(tags.length);
            for (int i = 0; i < this.getTagCount(); i++) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    dos.writeLong(tags[i][j]);
                }
                dos.writeByte(variableTagLengthInLong[i]);
                dos.writeShort(tagLength[i]);
                dos.writeBoolean(ifContig[i]);
                dos.writeInt(pairIndex[i]);
                dos.writeByte(mappingNum[i]);
                for (int j = 0; j < mappingNum[i]; j++) {
                    dos.writeByte(chr[i][j]);
                    dos.writeByte(strand[i][j]);
                    dos.writeInt(startPos[i][j]);
                    dos.writeShort(score[i][j]);
                    dos.writeByte(divergence[i][j]);
                }
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("File written to " + outfileS);
    }
    
    /**
     * read from a file
     * @param infileS 
     */
    public void readBinaryFile (String infileS) {
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infileS), 65536));
            tagLengthInLong = dis.readInt();
            int tagNum = dis.readInt();
            this.iniMatrix(tagNum);    
            for (int i = 0; i < tagNum; i++) {
                tags[i] = new long[tagLengthInLong];
                for (int j = 0; j < tagLengthInLong; j++) {
                    tags[i][j] = dis.readLong();
                }
                variableTagLengthInLong[i] = dis.readByte();
                tagLength[i] = dis.readShort();
                ifContig[i] = dis.readBoolean();
                pairIndex[i] = dis.readInt();
                mappingNum[i] = dis.readByte();
                if (mappingNum[i] == 0) continue;
                chr[i] = new byte[mappingNum[i]];
                strand[i] = new byte[mappingNum[i]];
                startPos[i] = new int[mappingNum[i]];
                score[i] = new short[mappingNum[i]];
                divergence[i] = new byte[mappingNum[i]];
                for (int j = 0; j < mappingNum[i]; j++) {
                    chr[i][j] = dis.readByte();
                    strand[i][j] = dis.readByte();
                    startPos[i][j] = dis.readInt();
                    score[i][j] = dis.readShort();
                    divergence[i][j] = dis.readByte();
                }
            }
            dis.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("File read from " + infileS);
    }
    
    private void importSamFile (String fastaFileS, String samFileS, int startIndex) {
        try {
            BufferedReader br = new BufferedReader (new FileReader(fastaFileS), 65536);
            for (int i = startIndex; i < startIndex+tags.length/2; i++) {
                if (br.readLine().split("_")[2].startsWith("c")) ifContig[i] = true;
                String seq = br.readLine();
                tagLength[i] = (short)seq.length();
                int left = tagLength[i] % BaseEncoder.chunkSize;
                if (left == 0) {
                    variableTagLengthInLong[i] = (byte)(tagLength[i]/BaseEncoder.chunkSize);
                }
                else {
                    variableTagLengthInLong[i] = (byte)(tagLength[i]/BaseEncoder.chunkSize+1);
                    StringBuilder sb = new StringBuilder();
                    for (int j = 0; j < BaseEncoder.chunkSize-left; j++) {
                        sb.append("A");
                    }
                    seq = seq + sb.toString();
                }
                tags[i] = BaseEncoder.getLongArrayFromSeq(seq);
                if (startIndex == 0) {
                    pairIndex[i] = i+tags.length/2;
                }
                else {
                    pairIndex[i] = i-tags.length/2;
                }     
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        try {
            BufferedReader br;
            if (samFileS.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(samFileS)))));
            } else {
                br = new BufferedReader(new FileReader(new File(samFileS)), 65536);
            }
            String temp;
            ArrayList<String> recordList = new ArrayList();
            while (!(temp = br.readLine()).startsWith("@PG")) {}
            temp = br.readLine();
            String[] tem = temp.substring(0, 30).split("\t");
            String currendId = tem[0];
            recordList.add(temp);
            while ((temp = br.readLine()) != null) {
                tem = temp.substring(0, 30).split("\t");
                if (tem[0].equals(currendId)) {
                    recordList.add(temp);
                }
                else {
                    importSAMRecord(recordList.toArray(new String[recordList.size()]), startIndex);
                    recordList = new ArrayList();
                    recordList.add(temp);
                    currendId = tem[0];
                }
            }
            importSAMRecord(recordList.toArray(new String[recordList.size()]), startIndex);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("SAM file was imported from " + samFileS);
    }
    
    private void importSAMRecord (String[] record, int startIndex) {
        int index = Integer.parseInt(record[0].split("_")[0])+startIndex;
        mappingNum[index] = (byte)record.length;
        chr[index] = new byte[record.length];
        strand[index] = new byte[record.length];
        startPos[index] = new int[record.length];
        score[index] = new short[record.length];
        divergence[index] = new byte[record.length];
        for (int i = 0; i < record.length; i++) {
            String[] temp = record[i].split("\\s");
            int orientiation=Integer.parseInt(temp[1]);
            int chromosome = Integer.MIN_VALUE;
            byte stran = Byte.MIN_VALUE;
            int sPos = Integer.MIN_VALUE;
            int ePos = Integer.MIN_VALUE;
            short mappingScore = Short.MIN_VALUE;
            byte diver = Byte.MIN_VALUE;
            if (orientiation == 4) {
                
                mappingNum[index] = 0;
                chr[index] = null;
                strand[index] = null;
                startPos[index] = null;
                score[index] = null;
                divergence[index] = null;
                continue;

            }
            else if (orientiation == 16 || orientiation == 272) {
                chromosome = Integer.parseInt(temp[2]);
                stran = -1;
                int[] alignSpan = SAMUtils.adjustCoordinates(temp[5], Integer.parseInt(temp[3]));
                sPos = alignSpan[1];
                ePos = alignSpan[0];
                mappingScore = Short.parseShort(temp[11].split(":")[2]);
                if (temp[17].startsWith("NM")) {
                    diver = Byte.parseByte(temp[17].split(":")[2]);
                }
                else {
                    diver = Byte.parseByte(temp[16].split(":")[2]);
                }
            }
            else {
                chromosome = Integer.parseInt(temp[2]);
                stran = 1;
                int[] alignSpan = SAMUtils.adjustCoordinates(temp[5], Integer.parseInt(temp[3]));
                sPos = alignSpan[0];
                ePos = alignSpan[1];
                mappingScore = Short.parseShort(temp[11].split(":")[2]);
                if (temp[17].startsWith("NM")) {
                    diver = Byte.parseByte(temp[17].split(":")[2]);
                }
                else {
                    diver = Byte.parseByte(temp[16].split(":")[2]);
                }
            }
            chr[index][i] = (byte)chromosome;
            strand[index][i] = stran;
            startPos[index][i] = sPos;
            score[index][i] = mappingScore;
            divergence[index][i] = diver;
        }
    }
    
    private void iniMatrix (int tagNum) {
        tags = new long[tagNum][];
        variableTagLengthInLong = new byte[tagNum];
        tagLength = new short[tagNum];
        ifContig = new boolean[tagNum];
        pairIndex = new int[tagNum];
        mappingNum = new byte[tagNum];
        chr = new byte[tagNum][];
        strand = new byte[tagNum][];
        startPos = new int[tagNum][];
        score = new short[tagNum][];
        divergence = new byte[tagNum][];
    }

    @Override
    public int getTagSizeInLong() {
        return tagLengthInLong;
    }

    public int getVariableTagSizeInLong (int index) {
        return this.variableTagLengthInLong[index];
    }
    
    @Override
    public String getNullTag() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getTagLength(int index) {
        return this.tagLength[index];
    }

    @Override
    public long[] getTag(int index) {
        return tags[index];
    }
    
    /**
     * Return the index of the tag which is the other end of PE
     * @param index
     * @return 
     */
    public int getPairIndex (int index) {
        return pairIndex[index];
    }
    
    /**
     * Return the number of mapping (multiple alignments)
     * @param index
     * @return 
     */
    public int getMappingNum (int index) {
        return mappingNum[index];
    }
    
    /**
     * Return chromosome of a mapping
     * @param tagIndex
     * @param mappingIndex
     * @return 
     */
    public byte getChr (int tagIndex, int mappingIndex) {
        return chr[tagIndex][mappingIndex];
    }
    
    /**
     * Return strand of a mapping
     * @param tagIndex
     * @param mappingIndex
     * @return 
     */
    public byte getStrand (int tagIndex, int mappingIndex) {
        return strand[tagIndex][mappingIndex];
    }
    
    /**
     * Return start position of a mapping
     * @param tagIndex
     * @param mappingIndex
     * @return 
     */
    public int getStartPos (int tagIndex, int mappingIndex) {
        return startPos[tagIndex][mappingIndex];
    } 
    
    /**
     * Return score of a mapping
     * @param tagIndex
     * @param mappingIndex
     * @return 
     */
    public short getScore (int tagIndex, int mappingIndex) {
        return score[tagIndex][mappingIndex];
    } 
    
    /**
     * Return divergence of a mapping
     * @param tagIndex
     * @param mappingIndex
     * @return 
     */
    public byte getDivergence (int tagIndex, int mappingIndex) {
        return divergence[tagIndex][mappingIndex];
    }

    /**
     * This may return one of multiple identical 64 bp tags, which are essentially differnt in full-length PE tags
     * @param read
     * @return 
     */
    @Override
    public int getTagIndex(long[] read) {
        //code inspired by COLT lower bound function
        int first = 0;
        int len = tags.length - first;
        int comp = 0;
        while (len > 0) {
            int half = len / 2;
            int middle = first + half;
            if ((comp = compareTags(middle, read)) < 0) {
                first = middle + 1;
                len -= half + 1;
            } else {
                len = half;
            }
        }
        if ((first < tags.length) && (compareTags(first, read) == 0)) {
            return first;
        }
        return -(first + 1);
    }

    /**
     * 
     * @param index
     * @param t2
     * @return 
     */
    private int compareTags(int index, long[] t2) {
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tags[index][i] < t2[i]) {
                return -1;
            }
            if (tags[index][i] > t2[i]) {
                return 1;
            }
        }
        return 0;
    }
    
    /**
     * Return the index of 64 bp tag with alignment of longest PE tag, when there are multiple identical 64 bp tags. When the tag is not found, return -1
     * @param read
     * @return The index of longest PE
     */
    public int getTagIndexWithLongestSeq (long[] read) {
        int[] range = this.getTagIndexSet(read);
        if (range == null) return -1;
        int index = -1, length = -1;
        for (int i = range[0]; i < range[1]; i++) {
            if (this.getTagLength(i) > length) {
                index = i;
                length = this.getTagLength(i);
            }
        }
        return index;
    }
    
    /**
     * Return the range of identical 64 bp tags, [0,3] means {0,1,2} are identical. When the read doesn't exist, return null
     * @param read
     * @return 
     */
    @Override
    public int[] getTagIndexSet(long[] read) {
        int index = this.getTagIndex(read);
        if (index < 0) return null;
        while (index > 0 && this.compareTags(index-1, read) == 0) index--;
        int[] range = new int[2];
        range[0] = index;
        while (index < this.tags.length-1 && this.compareTags(index+1, read) == 0) index++;
        range[1] = index+1;
        return range;
    }

    @Override
    public boolean areTagsUnique() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getTagCount() {
        return tags.length;
    }

    @Override
    public void swap(int index1, int index2) {
        long[] tmLong = tags[index1];
        tags[index1] = tags[index2];
        tags[index2] = tmLong;
        byte tmByte = variableTagLengthInLong[index1];
        variableTagLengthInLong[index1] = variableTagLengthInLong[index2];
        variableTagLengthInLong[index2] = tmByte;
        short tmShort = tagLength[index1];
        tagLength[index1] = tagLength[index2];
        tagLength[index2] = tmShort;
        boolean tmBoo = ifContig[index1];
        ifContig[index1] = ifContig[index2];
        ifContig[index2] = tmBoo;
        int tmInt1 = pairIndex[index1];
        int tmInt2 = pairIndex[index2];
        pairIndex[pairIndex[index1]] = index2;
        pairIndex[pairIndex[index2]] = index1;
        pairIndex[index1] = tmInt2;
        pairIndex[index2] = tmInt1;
        tmByte = mappingNum[index1];
        mappingNum[index1] = mappingNum[index2];
        mappingNum[index2] = tmByte;
        byte[] tmBytes = chr[index1];
        chr[index1] = chr[index2];
        chr[index2] = tmBytes;
        tmBytes = strand[index1];
        strand[index1] = strand[index2];
        strand[index2] = tmBytes;
        int[] tmInts = startPos[index1];
        startPos[index1] = startPos[index2];
        startPos[index2] = tmInts;
        short[] tmShorts = score[index1];
        score[index1] = score[index2];
        score[index2] = tmShorts;
        tmBytes = divergence[index1];
        divergence[index1] = divergence[index2];
        divergence[index2] = tmBytes;
        
    }

    @Override
    public int compare(int index1, int index2) {
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tags[index1][i] < tags[index2][i]) {
                return -1;
            }
            if (tags[index1][i] > tags[index2][i]) {
                return 1;
            }
        }
        return 0;
    }
    
    /**
     * Sort the class based on 64 bp tags
     */
    public void sort() {
        System.out.println("Position index sort begin.");
        GenericSorting.quickSort(0, this.getTagCount(), this, this);
        System.out.println("Position index sort end.");
    }
}
