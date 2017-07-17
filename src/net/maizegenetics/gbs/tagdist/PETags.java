package net.maizegenetics.gbs.tagdist;

import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;

/**
 * Basic interface for holding sets of sequence tag (these are compressed into 2-bit codings
 * of tags).  The length of tag in good sequence is also tracked.
 * @author Fei Lu
 */
public interface PETags extends Swapper, IntComparator {

    /**
     * Returns the number of long use to represent the sequence
     * @return Tag Length in Long primitive data type
     */
    public int getTagSizeInLong();

    /**
     * Returns the length in bp of a forward tag
     * @param index
     * @return length of the forward tag
     */
    public short getTagFLength(int index);
    
    /**
     * Returns the length in bp of a backward tag
     * @param index
     * @return length of the backward tag
     */
    public short getTagBLength(int index);

    /**
     * Get the compressed forward tag sequence in a long array for a given index
     * @param index
     * @return compressed forward tag sequence in long array
     */
    public long[] getTagF(int index);
    
    /**
     * Get the compressed backward tag sequence in a long array for a given index
     * @param index
     * @return compressed backward tag sequence in long array
     */
    public long[] getTagB(int index);

    /**
     * Gets the index of a PE tag (including forward and backward) (the only one if a unique list).
     * If the read is not found then it return
     * a negative value indicating its insertion point.
     * @param tagF
     *        Long array of forward tag
     * @param tagB
     *        Long array of backward tag
     * @return index of the PE tag
     */
    public int getTagIndex(long[] tagF, long[] tagB);


    /**
     * This is the number of different PE tags in the list (NOT THE SUM OF THE COUNTS)
     * The index will vary from 0 to (ReadTotal-1)
     * This is the number of distinct tags if readUnique is true
     * @return total number of tags
     */
    public int getTagCount();
    
    /**
     * Return the contig of PE tag
     * @param index
     * @return Contig of a PE tag, if there is no contig, return null
     */
    public long[] getContig (int index);
    
    /**
     * Return total number of PE contigs
     * @return Total count of PE Contig
     */
    
    public int getContigCount ();
    
    /**
     * Return the length of a contig
     * @param index
     * @return PE contig length
     */
    public short getContigLength (int index);
    
    /**
     * Return the contig length in Long primitive data type
     * @param index
     * @return contig length in Long primitive data type
     */
    public byte getContigLengthInLong (int index);
    
    
}
