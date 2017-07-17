/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.List;
import net.maizegenetics.gbs.maps.TagMappingInfoV3;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMapV3;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteHDF5TagGroups;

/**
 * Utils for tag genetic mapping
 * @author Fei Lu
 */
public class TagMappingUtils {
    
    public TagMappingUtils () {
        
    }
    
    public void mkTBTTagBlockFile (String tbtH5FileS, String topmH5FileS, String blockFileS, String software) {
        long lastTimePoint = System.nanoTime();
        TagsByTaxaByteHDF5TagGroups tbt = new TagsByTaxaByteHDF5TagGroups (tbtH5FileS);
        System.out.println("Loading TBT HDF5 took " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds");
        System.out.println("TBT has " + tbt.getTagCount() + " tags and " + tbt.getTaxaCount() + " taxa\n");
        lastTimePoint = System.nanoTime();
        TagMappingInfoV3.Aligner alignerName = TagMappingInfoV3.Aligner.getAlignerFromName(software);
        if (alignerName == null) {
            System.out.println("Input software is not Bowtie2, BWA or Blast, not supporting other aligner for now.");
            System.out.println("Program stops.");
            System.exit(0);
        }
        int[] blockChr = new int[tbt.getTagCount()];
        int[] blockPos = new int[tbt.getTagCount()];
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmH5FileS);
        long[] t;
        int index;
        for (int i = 0; i < tbt.getTagCount(); i++) {
            t = tbt.getTag(i);
            index = topm.getTagIndex(t);
            if (index < 0) {
                blockChr[i] = Integer.MIN_VALUE;
                blockPos[i] = Integer.MIN_VALUE;
            }
            else {
                int[] chrPos = topm.getUniqueMappingOfAligner(index, alignerName);
                if (chrPos == null) {
                    blockChr[i] = Integer.MIN_VALUE;
                    blockPos[i] = Integer.MIN_VALUE;
                }
                else {
                    blockChr[i] = chrPos[0];
                    blockPos[i] = chrPos[1];
                }
            }
        }
        System.out.println("Loading blocking mapping information took " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds\n");
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(blockFileS), 65536));
            dos.writeInt(blockChr.length);
            for (int i = 0; i < blockChr.length; i++) {
                dos.writeInt(blockChr[i]);
                dos.writeInt(blockPos[i]);
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    public ArrayList<int[]> getTBTTagBlock (String blockFileS) {
        ArrayList<int[]>  chrPosList = new ArrayList();
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(blockFileS), 65536));
            int[] blockChr, blockPos;
            int tagNum = dis.readInt();
            blockChr = new int[tagNum];
            blockPos = new int[tagNum];
            for (int i = 0; i < tagNum; i++) {
                blockChr[i] = dis.readInt();
                blockPos[i] = dis.readInt();
            }
            dis.close();
            chrPosList.add(blockChr);
            chrPosList.add(blockPos);
            return chrPosList;
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
        return chrPosList;
    }
    
    private double getTimeSpanSecond (long lastTimePoint) {
        return (double)this.getTimeSpanNano(lastTimePoint)/1000000000;
    }
    
     private long getTimeSpanNano (long lastTimePoint) {
        return System.nanoTime()- lastTimePoint;
    }
}
