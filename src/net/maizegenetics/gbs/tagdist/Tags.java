package net.maizegenetics.gbs.tagdist;

import cern.colt.Swapper;
import cern.colt.function.IntComparator;

/**
 * Basic interface for holding sets of sequence tag (these are compressed into 2-bit codings
 * of tags).  The length of tag in good sequence is also tracked.
 * <p>
 * Tags are encoded in longs with a 2-bit encoding see the package description for details
 * {@link net.maizegenetics.gbs.tagdist}.
 * <p>
 * @author Ed Buckler
 */
public interface Tags extends Swapper, IntComparator {

    /**
     * Returns the number of longs use to represent the sequence
     * @return
     */
    public int getTagSizeInLong();

    /**
     * Returns a polyA string of length 32*getTagSizeInLong() 
     * @return
     */
    public String getNullTag();

    /**
     * Returns the length in bp of a particular tag
     * @param index
     * @return length
     */
    public int getTagLength(int index);

    /**
     * Get the compressed read sequence in a long array for a given index
     * @param index
     * @return compressed read sequence in long array
     */
    public long[] getTag(int index);

    /**
     * Gets the first index of a read (the only one if a unique list).
     * If the read is not found then it return
     * a negative value indicating its insertion point.
     * @param read as a compressed long array
     * @return index of the read in the array
     */
    public int getTagIndex(long[] read);

    /**
     * Gets the set indices of matching reads (the only one if a unique list).
     * If the read is not found then it returns a null array
     * indicating its insertion point.
     * @param read as a compressed long array
     * @return set of indices of the read in the array (or null)
     */
    public int[] getTagIndexSet(long[] read);

    /**
     * Reports whether this list of reads includes duplicates
     * True is there are no read duplicates, false otherwise.
     * Collapses read files are likely to be unique
     * While a virtual digest of a genome with contain some duplicates
     * @return whether tags are unique
     */
    public boolean areTagsUnique();

    /**
     * This is the number of different tags in the list (NOT THE SUM OF THE COUNTS)
     * The index will vary from 0 to (ReadTotal-1)
     * This is the number of distinct tags if readUnique is true
     * @return total number of tags
     */
    public int getTagCount();
}
