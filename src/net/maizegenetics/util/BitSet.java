/*
 * BitSet
 */
package net.maizegenetics.util;

/**
 *
 * @author terry
 */
public interface BitSet {

    /**
     * Returns capacity.
     * 
     * @return capacity
     */
    public long capacity();

    /**
     * Returns the current capacity of this set.  Included for
     * compatibility.  This is *not* equal to {@link #cardinality}
     */
    public long size();

    /**
     * Returns true if there are no set bits
     */
    public boolean isEmpty();

    /**
     * Expert: returns the long[] storing the bits
     */
    public long[] getBits();
    
    /**
     * Expert: returns the long[] storing the bits from start to end
     */
    public long[] getBits(int startWord, int endWord);
    
    /**
     * Expert: returns 64 bits at index.
     */
    public long getBits(int index);

    /**
     * Expert: sets a new long[] to use as the bit storage
     */
    public void setBits(long[] bits);
    
    /**
     * Expert: sets specified word with given bits.
     * 
     * @param wordNum word index
     * @param bits bits
     */
    public void setLong(int wordNum, long bits);

    /**
     * Expert: gets the number of longs in the array that are in use
     */
    public int getNumWords();

    /**
     * Expert: sets the number of longs in the array that are in use
     */
    public void setNumWords(int nWords);

    /**
     * Returns true or false for the specified bit index.
     */
    public boolean get(int index);

    /**
     * Returns true or false for the specified bit index.
     * The index should be less than the BitSet size
     */
    public boolean fastGet(int index);

    /**
     * Returns true or false for the specified bit index
     */
    public boolean get(long index);

    /**
     * Returns true or false for the specified bit index.
     * The index should be less than the BitSet size.
     */
    public boolean fastGet(long index);

    /**
     * Returns 1 if the bit is set, 0 if not.
     * The index should be less than the BitSet size
     */
    public int getBit(int index);

    /**
     * Sets a bit, expanding the set size if necessary
     */
    public void set(long index);

    /**
     * Sets the bit at the specified index.
     * The index should be less than the BitSet size.
     */
    public void fastSet(int index);

    /**
     * Sets the bit at the specified index.
     * The index should be less than the BitSet size.
     */
    public void fastSet(long index);

    /**
     * Sets a range of bits, expanding the set size if necessary
     *
     * @param startIndex lower index
     * @param endIndex one-past the last bit to set
     */
    public void set(long startIndex, long endIndex);

    /**
     * Clears a bit.
     * The index should be less than the BitSet size.
     */
    public void fastClear(int index);

    /**
     * Clears a bit.
     * The index should be less than the BitSet size.
     */
    public void fastClear(long index);

    /**
     * Clears a bit, allowing access beyond the current set size without changing the size.
     */
    public void clear(long index);

    /**
     * Clears a range of bits.  Clearing past the end does not change the size of the set.
     *
     * @param startIndex lower index
     * @param endIndex one-past the last bit to clear
     */
    public void clear(int startIndex, int endIndex);

    /**
     * Clears a range of bits.  Clearing past the end does not change the size of the set.
     *
     * @param startIndex lower index
     * @param endIndex one-past the last bit to clear
     */
    public void clear(long startIndex, long endIndex);

    /**
     * Sets a bit and returns the previous value.
     * The index should be less than the BitSet size.
     */
    public boolean getAndSet(int index);

    /**
     * Sets a bit and returns the previous value.
     * The index should be less than the BitSet size.
     */
    public boolean getAndSet(long index);

    /**
     * Flips a bit.
     * The index should be less than the BitSet size.
     */
    public void fastFlip(int index);

    /**
     * Flips a bit.
     * The index should be less than the BitSet size.
     */
    public void fastFlip(long index);

    /**
     * Flips a bit, expanding the set size if necessary
     */
    public void flip(long index);

    /**
     * Flips a bit and returns the resulting bit value.
     * The index should be less than the BitSet size.
     */
    public boolean flipAndGet(int index);

    /**
     * Flips a bit and returns the resulting bit value.
     * The index should be less than the BitSet size.
     */
    public boolean flipAndGet(long index);

    /**
     * Flips a range of bits, expanding the set size if necessary
     *
     * @param startIndex lower index
     * @param endIndex one-past the last bit to flip
     */
    public void flip(long startIndex, long endIndex);

    /**
     * @return the number of set bits
     */
    public long cardinality();
    
    /**
     * Return number of set bits up to and including
     * bit at given index.
     * 
     * @param index index
     * 
     * @return the number of set bits
     */
    public long cardinality(int index);

    /**
     * Returns the index of the first set bit starting at the index specified.
     * -1 is returned if there are no more set bits.
     */
    public int nextSetBit(int index);

    /**
     * Returns the index of the first set bit starting at the index specified.
     * -1 is returned if there are no more set bits.
     */
    public long nextSetBit(long index);

    /**
     * Returns the index of the previous set bit starting at the index specified.
     * If the bit at index is set, index is returned, otherwise the next lower numbered set bit is returned.
     * -1 is returned if there are no more set bits.
     */
    public int previousSetBit(int index);

    /**
     * Returns the index of the previous set bit starting at the index specified.
     * If the bit at index is set, index is returned, otherwise the next lower numbered set bit is returned.
     * -1 is returned if there are no more set bits.
     */
    public long previousSetBit(long index);

    /**
     * this = this AND other
     */
    public void intersect(BitSet other);

    /**
     * this = this OR other
     */
    public void union(BitSet other);

    /**
     * Remove all elements set in other. this = this AND_NOT other
     */
    public void remove(BitSet other);

    /**
     * XOR
     */
    public void xor(BitSet other);

    /**
     * AND
     *
     * @param other
     */
    public void and(BitSet other);

    /**
     * OR
     */
    public void or(BitSet other);

    /**
     * see {@link andNot}
     */
    public void andNot(BitSet other);

    /**
     * Returns true if the sets have any elements in common
     */
    public boolean intersects(BitSet other);

    /**
     * Expand the long[] with the size given as a number of words (64 bit longs).
     * getNumWords() is unchanged by this call.
     */
    public void ensureCapacityWords(int numWords);

    /**
     * Ensure that the long[] is big enough to hold numBits, expanding it if necessary.
     * getNumWords() is unchanged by this call.
     */
    public void ensureCapacity(long numBits);

    /**
     * Lowers numWords, the number of words in use,
     * by checking for trailing zero words.
     */
    public void trimTrailingZeros();

    /**
     * Returns index of the nth set bit.
     *
     * @param n nth set bit
     *
     * @return index
     */
    public int indexOfNthSetBit(int n);

    /**
     * Return indices of set bits.
     *
     * @return indices
     */
    public int[] getIndicesOfSetBits();
}
