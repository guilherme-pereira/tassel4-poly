/*
 * TagsByTaxa
 */
package net.maizegenetics.gbs.tagdist;

import net.maizegenetics.util.OpenBitSet;

import java.io.File;

/**
 * Tags by Taxa is the interface for tracking the distribution of tags across
 * taxa.  It can be thought of as a matrix of nTaxa by nTags.  Taxa names are
 * also maintained by this interface.
 *
 *
 * @author Ed Buckler, Gabriel Rodrigues Alves Margarido
 */
public interface TagsByTaxa extends Tags {

    /**
     * Types of file packing used by GBS data structures.  Bit compresses all information
     * to presence/absence, while Byte/Short/Int tradeoff memory size with depth of coverage.
     * Text is only used for file storage.
     */
    public static enum FilePacking {

        Bit, Byte, Short, Int, Text
    };

    /**
     * Method for incrementing the read depth for the specified tag and taxon.
     * @param tagIndex
     * @param taxaIndex
     * @param addValue value to add to the current value
     */
    void addReadsToTagTaxon(int tagIndex, int taxaIndex, int addValue);

    /**
     * Get the index of a taxon based on its name
     * @param taxon name of the taxon (sample)
     * @return  taxon index
     */
    int getIndexOfTaxaName(String taxon);

    /**
     * Returns the read count for a specified tag and taxon
     * @param tagIndex
     * @param taxaIndex
     * @return read depth
     */
    int getReadCountForTagTaxon(int tagIndex, int taxaIndex);

    /**
     * Returns all read counts for a specified tag index.
     * @param tagIndex  index of the tag
     * @return array of read depths in order of taxa
     */
    short[] getTaxaReadCountsForTag(int tagIndex);

    /**
     * TODO: Jeff annotate this
     */
    public void truncateTaxonNames();


    /**
     *
     * @param rowSetMethod
     */
    public void setMethodByRows(boolean rowSetMethod);

    /**
     * The presence/absence of the taxa in a bitSet format (1=present, 0=absent)
     * @param tagIndex the index of the tag
     * @return BitSet representing distribution of taxa
     */
    OpenBitSet getTaxaReadBitsForTag(int tagIndex);

    /**
     * Returns the number of taxa (samples)
     * @return
     */
    int getTaxaCount();

    /**
     * Get the number of taxa with a given tag 
     * It is count of the taxa with readCount>0
     * @param readIndex
     * @return number of taxa with a read
     */
    int getNumberOfTaxaWithTag(int readIndex);

    /**
     * Returns the name of a taxon
     * @param taxaIndex the index of the taxon
     * @return taxon name
     */
    String getTaxaName(int taxaIndex);

    /**
     * Returns an array of taxa names
     */
    String[] getTaxaNames();

    /**
     * Set the read depth for a specific tag and taxon.  Overwrite the depth that was previously there.
     * @param tagIndex index of the tag
     * @param taxaIndex  index of the taxon
     * @param value  read depth
     */
    void setReadCountForTagTaxon(int tagIndex, int taxaIndex, int value);

    /**
     * Write the TagsByTaxa object to a file in one of the FilePacking formats.  MinCount is useful for excluding the low
     * depth tags (many of these are sequencing errors).
     * @param outFile File for saving this object
     * @param fileType  Determine how to output the depth
     * @param minCount  Only tags greater than minCount are output
     */
    void writeDistFile(File outFile, FilePacking fileType, int minCount);

    /**Initializes the storage matrices for reading (unclear why in interface)*/
    void initMatrices(int taxaNum, int tagNum);

    /**
     * Add taxa to the TagsByTaxa matrix, they will all be set to distribution value of
     * zero.
     * @param addTaxaNames
     */
    void addTaxa(String[] addTaxaNames);

    /**In implementations that use a RandomAccessFile for storage, this clears the RAM buffer of any
    remaining data, writes it to the file on disk, and closes the file.*/
    public void getFileReadyForClosing();
}
