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
import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.util.BaseEncoder;

/**
 * Class that hold genetic position of tags, genetic positions are the position of markers associated with tags
 * @author Fei Lu
 */
public class TagsOnGeneticMap extends AbstractTags {
    protected int[] gChr;
    protected int[] gPos;
    
    /**
     * Construct TOGM from a file
     * @param infileS
     *        File name of input TagsOnGeneticMap file
     * @param format 
     *        FilePacking format
     */
    public TagsOnGeneticMap (String infileS, FilePacking format) {
        this.readDistFile(infileS, format);
    }
    
    /**
     * Return chromosome of genetic position
     * @param index
     * @return Chromosome of genetic position 
     */
    public int getGChr (int index) {
        return gChr[index];
    }
    
    /**
     * Return site of genetic position
     * @param index
     * @return Genetic position 
     */
    public int getGPos (int index) {
        return gPos[index];
    }
    
    /**
     * Read tagsOnGeneticMap file
     * @param infileS
     * @param format 
     */
    public void readDistFile (String infileS, FilePacking format) {
        System.out.println("Reading TOGM file to " + infileS);
        File infile = new File (infileS);
        switch (format) {
            case Text:
                readTextTOGMFile(infile);
                break;
            default:
                readBinaryTOGMFile(infile);
                break;
        }
        System.out.println("TOGM file read. Tatol: " + this.getTagCount() + " PETags");
    }
    
    /**
     * Read text TOGM file
     * @param infile 
     */
    private void readTextTOGMFile (File infile) {
        try {
            BufferedReader br = new BufferedReader (new FileReader(infile), 65536);
            tagLengthInLong = Integer.parseInt(br.readLine());
            this.iniMatrix(tagLengthInLong, Integer.parseInt(br.readLine()));
            br.readLine();
            for (int i = 0; i < this.getTagCount(); i++) {
                String[] temp = br.readLine().split("\t");
                long[] t = BaseEncoder.getLongArrayFromSeq(temp[0]);
                for (int j = 0; j < tagLengthInLong; j++) {
                    tags[j][i] = t[j];
                }
                tagLength[i] = Byte.parseByte(temp[1]);
                gChr[i] = Integer.parseInt(temp[2]);
                gPos[i] = Integer.parseInt(temp[3]);        
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Read binary TOGM file
     * @param infile 
     */
    private void readBinaryTOGMFile (File infile) {
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infile), 65536));
            tagLengthInLong = dis.readInt();
            this.iniMatrix(tagLengthInLong, dis.readInt());
            for (int i = 0; i < this.getTagCount(); i++) {
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    this.tags[j][i] = dis.readLong();
                }
                this.tagLength[i] = dis.readByte();
                this.gChr[i] = dis.readInt();
                this.gPos[i] = dis.readInt();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Initialize the matrix of TOGM
     * @param tagLengthInLong
     *        Tag length in Long primitive data type
     * @param tagNum 
     *        Total tag number
     */
    protected void iniMatrix (int tagLengthInLong, int tagNum) {
        tags = new long[tagLengthInLong][tagNum];
        tagLength = new byte[tagNum];
        gChr = new int[tagNum];
        gPos = new int[tagNum];
    }
    
    /**
     * Write TagsOnGeneticMap file
     * @param outfileS
     *        File name of output file
     * @param format 
     *        FilePacking format
     */
    public void writeDistFile (String outfileS, FilePacking format) {
        System.out.println("Writing TOGM file to " + outfileS);
        switch (format) {
            case Text:
                writeTextTOGMFile(outfileS);
                break;
            default:
                writeBinaryTOGMFile(outfileS);
                break;
        }
        System.out.println("TOGM file written");
    }
    
    /**
     * Write text TOGM file
     * @param outfileS 
     */
    private void writeTextTOGMFile (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write(String.valueOf(this.tagLengthInLong));
            bw.newLine();
            bw.write(String.valueOf(this.getTagCount()));
            bw.newLine();
            bw.write("Tag\tTagLength\tGChr\tGPos");
            bw.newLine();
            long[] temp = new long[this.tagLengthInLong];
            for (int i = 0; i < this.getTagCount(); i++) {
				for (int j = 0; j < temp.length; j++) {
					temp[j] = tags[j][i];
				}
                bw.write(BaseEncoder.getSequenceFromLong(temp)+"\t"+String.valueOf(this.getTagLength(i))+"\t");
                bw.write(String.valueOf(this.gChr[i])+"\t"+String.valueOf(this.gPos[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Write binary TOGM file
     * @param outfileS 
     */
    private void writeBinaryTOGMFile (String outfileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
            dos.writeInt(tagLengthInLong);
            dos.write(this.getTagCount());
            for (int i = 0; i < this.getTagCount(); i++) {
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    dos.writeLong(this.tags[j][i]);
                }
                dos.writeByte(this.getTagLength(i));
                dos.writeInt(this.getGChr(i));
                dos.writeInt(this.getGPos(i));
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }  
}
