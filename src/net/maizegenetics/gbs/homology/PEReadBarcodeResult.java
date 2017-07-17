/*
 * ReadBarcodeResult
 */
package net.maizegenetics.gbs.homology;

import net.maizegenetics.gbs.util.BaseEncoder;

/**
 * Container class for returning the results of parsed barcoded sequencing read.
 * 
 * @author Fei Lu
 */
public class PEReadBarcodeResult {
    /**Long array of forward sequence*/
    long[] readF;
    /**Length of forward sequence*/
    short lengthF;
    /**Long array of backward sequence*/
    long[] readB;
    /**Length of backward sequence*/
    short lengthB;
    /**Taxon name of the sequence*/
    String taxonName;

    /**
     * Construct object from forward and backward sequence
     * @param rF
     * @param rB 
     */
    public PEReadBarcodeResult(ShortReadBarcodeResult rF, ShortReadBarcodeResult rB) {
        this.taxonName = rF.taxonName;
        this.readF = rF.read;
        this.lengthF = rF.length;
        this.readB = rB.read;
        this.lengthB = rB.length;
        this.orderReadFReadB();
    }

    /**
     * Return length of forward sequence 
     * @return length of forward sequence 
     */
    public short getLengthF() {
        return lengthF;
    }

    /**
     * Return length of backward sequence 
     * @return length of backward sequence 
     */
    public short getLengthB() {
        return lengthB;
    }
    
    /**
     * Return forward sequence in Long array 
     * @return forward sequence in Long array
     */
    public long[] getReadF() {
        return readF;
    }
    
    /**
     * Return backward sequence in Long array 
     * @return backward sequence in Long array
     */
    public long[] getReadB() {
        return readB;
    }

    /**
     * Return taxon name of the sequence
     * @return taxon name of the sequence
     */
    public String getTaxonName() {
        return taxonName;
    }
    
    /**
     * Return tag length in Long primitive data type
     * @return tag length in Long primitive data type
     */
    public int getTagLengthInLong () {
        return readF.length;
    }
    
    /**
     * Order the forward sequence and backward sequence
     * @return boolean value of if the forward and backward sequence is switched 
     */
    private boolean orderReadFReadB () {
        if (this.compareReadFReadB() == 1) {
            this.switchReadFReadB();
            return true;
        }
        return false;
    }
    
    /**
     * Switch the forward sequence and backward sequence
     */
    private void switchReadFReadB () {
        for (int i = 0; i < readF.length; i++) {
            long temp = readF[i];
            readF[i] = readB[i];
            readB[i] = temp;  
        }
        short tem = lengthF;
        lengthF = lengthB;
        lengthB = tem;
    }
    
    /**
     * Compare the forward sequence and backward sequence
     * @return result of comparison -1(<), 1(>), 0(=) 
     */
    private int compareReadFReadB () {
        int tagLengthInLong = readF.length;
        for (int i = 0; i < tagLengthInLong; i++) {
            if (readF[i] < readB[i]) return -1;
            if (readF[i] > readB[i]) return 1;
        }
        return 0;
    }
}
