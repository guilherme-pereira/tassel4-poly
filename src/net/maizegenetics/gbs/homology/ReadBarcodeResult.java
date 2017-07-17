/*
 * ReadBarcodeResult
 */
package net.maizegenetics.gbs.homology;

import net.maizegenetics.gbs.util.BaseEncoder;

/**
 * Container class for returning the results of parsed barcoded sequencing read.
 * <p>
 * {@link #unprocessedSequence} is the original sequence with the barcode still attached
 * <p>
 * {@link #processedSequence} is the sequence with the barcode removed and cut down to the given {@link length}.
 * <p>
 * {@link #paddedSequence} is the {@link #processedSequence} padded with polyA to length if the 
 * {@link #processedSequence} was shorter than {@link length}.
 * 
 * @author Ed Buckler
 */
public class ReadBarcodeResult {
    /**Original sequence from sequencer*/
    public String unprocessedSequence = null;
    /**Processed sequence with barcode removed*/
    public String processedSequence = null;
    /**Processed sequence padded with polyA*/
    public String paddedSequence = null;
    /**length of the processed sequence*/
    byte length;
    /**Sequence encoded in 2-bit long array*/
    long[] read;
    /**Taxon name implied by the barcode sequence*/
    String taxonName;

    //TODO this instantiation should also include the orginal unprocessedSequence, processedSequence, and paddedSequence - the the object encode it
    public ReadBarcodeResult(long[] read, byte length, String taxon) {
        this.read = read;
        this.length = length;
        this.taxonName = taxon;
    }

    public ReadBarcodeResult(String sequence) {
        unprocessedSequence = sequence;
    }

    @Override
    public String toString() {
        return BaseEncoder.getSequenceFromLong(read) + ":" + (int) length + ":" + taxonName;
    }

    public byte getLength() {
        return length;
    }

    public long[] getRead() {
        return read;
    }

    /**Return taxon name implied by the barcode sequence*/
    public String getTaxonName() {
        return taxonName;
    }
}
