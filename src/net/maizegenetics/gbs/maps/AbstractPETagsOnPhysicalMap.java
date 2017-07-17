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
import net.maizegenetics.gbs.tagdist.AbstractPETags;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.util.BaseEncoder;

/**
 * Basic method of working with alignment of PE tags
 * @author Fei Lu
 * @deprecated 
 */
public abstract class AbstractPETagsOnPhysicalMap extends AbstractPETags {
    /***/
    protected int[] chrF;//Integer.MIN_VALUE means the tag is not physically mapped
    protected int[] chrB;
    protected int[] chrContig;
    protected int[] posStartF;//Integer.MIN_VALUE means the tag is not physically mapped
    protected int[] posStartB;
    protected int[] posStartContig;
    protected int[] posEndF;
    protected int[] posEndB;
    protected int[] posEndContig;
    protected byte[] strandF; //Byte.MIN_VALUE means the tag is not physically mapped
    protected byte[] strandB;
    protected byte[] strandContig;
    
    public enum PETagType {
        Forward, Backward, Contig
    }
    
    @Override
    protected void iniMatrix (int tagLengthInLong, int tagNum) {
        super.iniMatrix(tagLengthInLong, tagNum);
        chrF = new int[tagNum];
        chrB = new int[tagNum];
        chrContig = new int[tagNum];
        posStartF = new int[tagNum];
        posStartB = new int[tagNum];
        posStartContig = new int[tagNum];
        posEndF = new int[tagNum];
        posEndB = new int[tagNum];
        posEndContig = new int[tagNum];
        strandF = new byte[tagNum];
        strandB = new byte[tagNum];
        strandContig = new byte[tagNum];
    }
    
    public int getChrF (int index) {
        return chrF[index];
    }
    
    public int getChrB (int index) {
        return chrB[index];
    }
    
    public int getChrContig (int index) {
        return chrContig[index];
    }
    
    public int getPosStartF (int index) {
        return posStartF[index];
    }
    
    public int getPosStartB (int index) {
        return posStartB[index];
    }
    
    public int getPosStartContig (int index) {
        return posStartContig[index];
    }
    
    public int getPosEndF (int index) {
        return posEndF[index];
    }
    
    public int getPosEndB (int index) {
        return posEndB[index];
    }
    
    public int getPosEndContig (int index) {
        return posEndContig[index];
    }
    
    public byte getStrandF (int index) {
        return strandF[index];
    }
    
    public byte getStrandB (int index) {
        return strandB[index];
    }
    
    public byte getStrandContig (int index) {
        return strandContig[index];
    }
}
