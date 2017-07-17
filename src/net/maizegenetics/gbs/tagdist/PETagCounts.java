/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.tagdist;

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
import java.util.TreeSet;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.util.BaseEncoder;

/**
 * Hold PE tags and their counts, also including contigs
 * @author Fei Lu
 */
public class PETagCounts extends AbstractPETags {
    /** Read counts of PE tags*/
    int[] readCount;
    
    /**
     * Construct PETagCounts from a file
     * @param inFile
     *        Filename of PETagCounts
     * @param format 
     *        FilePacking format
     */
    public PETagCounts (String inFile, FilePacking format) {
        this.readDistFile(inFile, format);
    }
    
    /**
     * Initialize PETagCounts with empty matrix
     * @param tagLengthInLong
     *        Tag length In Long primitive data type
     * @param tagNum 
     *        Tag number
     */
    public PETagCounts (int tagLengthInLong, int tagNum) {
        this.iniMatrix(tagLengthInLong, tagNum);
    }
    
    /**
     * Collapse the PETagCounts and return a new collapsed PETagCounts object
     * @return Collapsed PETagCounts object 
     */
    public PETagCounts getCollapsedPETagCounts () {
        int tagNum = this.getTagCount() - this.collapseCounts();
        PETagCounts petc = new PETagCounts (this.getTagSizeInLong(), tagNum);
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (this.getReadCount(i) == 0) continue;
            for (int j = 0; j < this.getTagSizeInLong(); j++) {
                petc.tagsF[cnt][j] = this.tagsF[i][j];
                petc.tagsB[cnt][j] = this.tagsB[i][j];
            }
            petc.tagFLength[cnt] = this.tagFLength[i];
            petc.tagBLength[cnt] = this.tagBLength[i];
            petc.contigLengthInLong[cnt] = this.contigLengthInLong[i];
            petc.contig[cnt] = this.contig[i];
            petc.contigLength[cnt] = this.contigLength[i];
            petc.readCount[cnt] = this.readCount[i];
            cnt++;
        }
        return petc;
    }
    
    /**
     * Collapse PETagCounts, the tag count of collapsed tags is set to 0
     * @return total read count of collapsed PE tags  
     */
    public int collapseCounts () {
        this.sort();
        int collapsedRows = 0;
        for (int i = 1; i < this.getTagCount(); i++) {
            if (this.compare(i, i-1) == 0) {
                readCount[i] += readCount[i - 1];
                readCount[i - 1] = 0;
                collapsedRows++;
            }
        }
        System.out.println("Rows collapsed:" + collapsedRows);
        return collapsedRows;
    }
    
    /**
     * Merge two PETagCounts objects
     * @param another
     *        Another PETagCounts object
     * @param ifCollapsed
     *        Boolean value of another PETagCounts object (If it is collapsed). Both objects should be collapsed first
     * @return 
     */
    public PETagCounts getMergedPETagCounts (PETagCounts another, boolean ifCollapsed) {
        if (!ifCollapsed) another = another.getCollapsedPETagCounts();
        PETagCounts petc = new PETagCounts(this.tagLengthInLong, this.getTagCount()+another.getTagCount());
        int cnt = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            for (int j = 0; j < this.getTagSizeInLong(); j++) {
                petc.tagsF[cnt][j] = this.tagsF[i][j];
                petc.tagsB[cnt][j] = this.tagsB[i][j];
            }
            petc.tagFLength[cnt] = this.tagFLength[i];
            petc.tagBLength[cnt] = this.tagBLength[i];
            petc.contigLengthInLong[cnt] = this.contigLengthInLong[i];
            petc.contig[cnt] = this.contig[i];
            petc.contigLength[cnt] = this.contigLength[i];
            petc.readCount[cnt] = this.readCount[i];
            cnt++;
        }
        for (int i = 0; i < another.getTagCount(); i++) {
            for (int j = 0; j < another.getTagSizeInLong(); j++) {
                petc.tagsF[cnt][j] = another.tagsF[i][j];
                petc.tagsB[cnt][j] = another.tagsB[i][j];
            }
            petc.tagFLength[cnt] = another.tagFLength[i];
            petc.tagBLength[cnt] = another.tagBLength[i];
            petc.contigLengthInLong[cnt] = another.contigLengthInLong[i];
            petc.contig[cnt] = another.contig[i];
            petc.contigLength[cnt] = another.contigLength[i];
            petc.readCount[cnt] = another.readCount[i];
            cnt++;
        }
        return petc.getCollapsedPETagCounts();
    }
    
    /**
     * Return total read count of PE tags
     * @return Total read count of this PETagCounts object 
     */
    public int getTotalReadCount () {
        int sum = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            sum += this.getReadCount(i);
        }
        return sum;
    }
    
    /**
     * Return read count of a PE tag
     * @param index
     * @return Read count of a PE tag 
     */
    public int getReadCount (int index) {
        return readCount[index];
    }
    
    /**
     * Return total number of tags with a count greater than minimum count
     * @param minCount
     *        Minimum count of PE tag
     * @return Total number of tags with a count greater than minCount
     */
    private int getTagNumWithMincount (int minCount) {
        int num = 0;
        for (int i = 0; i < this.getTagCount(); i++) {
            if (readCount[i] >= minCount) num++;
        }
        return num;
    }
    
    /**
     * Read PETagCounts file
     * @param infileS
     *        File name of PETagCounts file
     * @param format 
     *        FilePacking format
     */
    public void readDistFile (String infileS, FilePacking format) {
        System.out.println("Reading PETagCounts file to " + infileS);
        File infile = new File (infileS);
        switch (format) {
            case Text:
                readTextPETagCountsFile(infile);
                break;
            default:
                readBinaryPETagCountsFile(infile);
                break;
        }
        System.out.println("PETagCounts file read. Tatol: " + this.getTagCount() + " PETags");
    }
    
    /**
     * Read binary PETagCounts file
     * @param infile 
     */
    private void readBinaryPETagCountsFile (File infile) {
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infile), 65536));
            tagLengthInLong = dis.readInt();
            int tagNum = dis.readInt();
            if (tagNum == -1) {
                int lineSize = (tagLengthInLong*8+2) * 2 + 1 + 2 + 4;
                tagNum = (int)((infile.length()-8)/lineSize);
            }
            this.iniMatrix(tagLengthInLong, tagNum);
            for (int i = 0; i < tagNum; i++) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    tagsF[i][j] = dis.readLong();
                }
                tagFLength[i] = dis.readShort();
                for (int j = 0; j < tagLengthInLong; j++) {
                    tagsB[i][j] = dis.readLong();
                }
                tagBLength[i] = dis.readShort();
                contigLengthInLong[i] = dis.readByte();
                contig[i] = new long[contigLengthInLong[i]];
                for (int j = 0; j < contig[i].length; j++) {
                    contig[i][j] = dis.readLong();
                }   
                contigLength[i] = dis.readShort();
                readCount[i] = dis.readInt();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Read text PETagCounts file
     * @param infile 
     *        File name of input file
     */
    private void readTextPETagCountsFile (File infile) {
        try {
            BufferedReader br = new BufferedReader (new FileReader(infile), 65536);
            tagLengthInLong = Integer.valueOf(br.readLine());
            int tagNum = Integer.valueOf(br.readLine()); 
            this.iniMatrix(tagLengthInLong, tagNum);
            for (int i = 0; i < tagNum; i++) {
                String[] temp = br.readLine().split("\t");
                long[] t = BaseEncoder.getLongArrayFromSeq(temp[0]);
                for (int j = 0; j < tagLengthInLong; j++) {
                    tagsF[i][j] = t[j];
                }
                tagFLength[i] = Short.valueOf(temp[1]);
                t = BaseEncoder.getLongArrayFromSeq(temp[2]);
                for (int j = 0; j < tagLengthInLong; j++) {
                    tagsB[i][j] = t[j];
                }
                tagBLength[i] = Short.valueOf(temp[3]);
                readCount[i] = Integer.valueOf(temp[4]);
                contigLengthInLong[i] = Byte.valueOf(temp[5]);
                contigLength[i] = Short.valueOf(temp[6]);
                this.contig[i] = new long[contigLengthInLong[i]];
                if (contigLengthInLong[i] != 0) {
                    t = BaseEncoder.getLongArrayFromSeq(temp[7]);
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
    
    /**
     * Write PETagCounts file
     * @param outfileS
     *        File name of output PETagCounts file
     * @param format
     *        FilePakcing format
     * @param minCount 
     *        Minimum count of PE Tag in the output
     */
    public void writeDistFile (String outfileS, FilePacking format, int minCount) {
        System.out.println("Writing PETagCounts file to " + outfileS);
        int outTagNum = this.getTagNumWithMincount(minCount);
        switch (format) {
            case Text:
                writeTextPETagCountsFile(outfileS, outTagNum, minCount);
                break;
            default:
                writeBinaryPETagCountsFile(outfileS, outTagNum, minCount);
                break;
        }
        System.out.println("PETagCounts file written");
    }
    
    /**
     * Write binary PETagCounts file
     * @param outfileS
     *        File name of output PETagCounts file
     * @param outTagNum
     *        Total number of PE tags in the output
     * @param minCount 
     *        Minimum count of PE Tag in the output
     */
    private void writeBinaryPETagCountsFile (String outfileS, int outTagNum, int minCount) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
            dos.writeInt(tagLengthInLong);
            dos.writeInt(outTagNum);
            for (int i = 0; i < this.getTagCount(); i++) {
                if (readCount[i] < minCount) continue;
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    dos.writeLong(this.tagsF[i][j]);
                }
                dos.writeShort(this.tagFLength[i]);
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    dos.writeLong(this.tagsB[i][j]);
                }
                dos.writeShort(this.tagBLength[i]);
                dos.writeByte(this.contigLengthInLong[i]);
                for (int j = 0; j < this.contigLengthInLong[i]; j++) {
                    dos.writeLong(contig[i][j]);
                }
                dos.writeShort(this.contigLength[i]);
                dos.writeInt(this.readCount[i]);
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
     * Write text PETagCounts file
     * @param outfileS
     *        File name of output PETagCounts file
     * @param outTagNum
     *        Total number of PE tags in the output
     * @param minCount 
     *        Minimum count of PE Tag in the output
     */
    private void writeTextPETagCountsFile (String outfileS, int outTagNum, int minCount) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write(String.valueOf(tagLengthInLong));
            bw.newLine();
            bw.write(String.valueOf(outTagNum));
            bw.newLine();
            for (int i = 0; i < this.getTagCount(); i++) {
                if (this.getReadCount(i) < minCount) continue;
                bw.write(BaseEncoder.getSequenceFromLong(this.getTagF(i))+"\t"+String.valueOf(this.getTagFLength(i))+"\t");
                bw.write(BaseEncoder.getSequenceFromLong(this.getTagB(i))+"\t"+String.valueOf(this.getTagBLength(i))+"\t");
                bw.write(String.valueOf(this.getReadCount(i))+"\t");
                bw.write(String.valueOf(this.contigLengthInLong[i])+"\t"+String.valueOf(this.contigLength[i])+"\t");
                bw.write(BaseEncoder.getSequenceFromLong(this.contig[i])+"\t");
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
    
    @Override
    protected void iniMatrix (int tagLengthInLong, int tagNum) {
        super.iniMatrix(tagLengthInLong, tagNum);
        readCount = new int[tagNum];
    }
    
    @Override
    public void swap(int index1, int index2) {
         super.swap(index1, index2);
         int temp = readCount[index1];
         readCount[index1] = readCount[index2];
         readCount[index2] = temp;
    }
}
