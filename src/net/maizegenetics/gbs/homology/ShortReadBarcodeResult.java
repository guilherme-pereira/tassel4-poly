/*
 * ReadBarcodeResult
 */
package net.maizegenetics.gbs.homology;

import net.maizegenetics.gbs.util.BaseEncoder;

/**
 * Container class for returning the results of parsed barcoded sequencing read.
 * The length of read is in short. Max length is 32767 bp.
 * 
 * @deprecated Fei should remove this class and use ReadBarcodeResult with length replaced
 * with short.
 * @author Fei Lu
 */
@Deprecated
public class ShortReadBarcodeResult {

    public String unprocessedSequence = null;
    public String processedSequence = null;
    public String paddedSequence = null;
    short length;
    long[] read;
    String taxonName;

    //TODO this instantiation should also include the orginal unprocessedSequence, processedSequence, and paddedSequence - the the object encode it
    public ShortReadBarcodeResult(long[] read, short length, String taxon) {
        this.read = read;
        this.length = length;
        this.taxonName = taxon;
    }
    
    /**
     * Construct object with genomic sequence
     * @param sequence 
     */
    public ShortReadBarcodeResult(String sequence) {
        unprocessedSequence = sequence;
    }

    @Override
    public String toString() {
        return BaseEncoder.getSequenceFromLong(read) + ":" + (int) length + ":" + taxonName;
    }

    /**
     * Return the length of the genomic sequence
     * @return Length of the genomic sequence
     */
    public short getLength() {
        return length;
    }

    public long[] getRead() {
        return read;
    }

    public String getTaxonName() {
        return taxonName;
    }
}
