/*
 * Reads
 */
package net.maizegenetics.gbs.util;

/**
 * Basic methods for working with nextgen reads
 * @author edbuckler
 */
public interface Reads {

    /**
     * Get the compressed read sequence in a long array for a given index
     * @param index
     * @return compressed read sequence in long array
     */
    public long[] getRead(int index);

    /**
     * Gets the count of a particular read.  If index does not exist it return -1
     * @param index
     * @return read count
     */
    public int getReadCount(int index);

    /**
     * Gets the first index of a read (the only one if a unique list).
     * If the read is not found then it return
     * a negative value indicating its insertion point.
     * @param read as a compressed long array
     * @return the index of the read in the array
     */
    public int getReadIndex(long[] read);

    /**
     * Gets the set indices of matching reads (the only one if a unique list).
     * If the read is not found then it returns a null array
     * indicating its insertion point.
     * @param read as a compressed long array
     * @return set of indices of the read in the array (or null)
     */
    public int[] getReadIndexSet(long[] read);

    /**
     * Reports whether this list of reads includes duplicates
     * True is there are no read duplicates, false otherwise.
     * Collapses read files are likely to be unique
     * While a virtual digest of a genome with contain some duplicates
     * @return
     */
    public boolean areReadsUnique();

    /**
     * This is the number of different reads in the list (NOT THE SUM OF THE COUNTS)
     * The index will vary from 0 to (ReadTotal-1)
     * This is the number of distinct reads if readUnique is true
     * @return total number of reads
     */
    public int getReadTotal();
}
