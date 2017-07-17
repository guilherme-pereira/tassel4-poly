/**
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements. See the NOTICE file distributed with this
 * work for additional information regarding copyright ownership. The ASF
 * licenses this file to You under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations under
 * the License.
 */
package net.maizegenetics.util;

import java.io.Serializable;
import java.util.Arrays;

//import org.apache.lucene.search.DocIdSetIterator;
/**
 * An "open" BitSet implementation that allows direct access to the array of
 * words storing the bits.
 * <p/>
 * Unlike java.util.bitset, the fact that bits are packed into an array of longs
 * is part of the interface. This allows efficient implementation of other
 * algorithms by someone other than the author. It also allows one to
 * efficiently implement alternate serialization or interchange formats.
 * <p/>
 * <
 * code>OpenBitSet</code> is faster than
 * <code>java.util.BitSet</code> in most operations and *much* faster at
 * calculating cardinality of sets and results of set operations. It can also
 * handle sets of larger cardinality (up to 64 * 2**32-1)
 * <p/>
 * The goals of
 * <code>OpenBitSet</code> are the fastest implementation possible, and maximum
 * code reuse. Extra safety and encapsulation may always be built on top, but if
 * that's built in, the cost can never be removed (and hence people re-implement
 * their own version in order to get better performance). If you want a "safe",
 * totally encapsulated (and slower and limited) BitSet class, use
 * <code>java.util.BitSet</code>.
 * <p/>
 * <h3>Performance Results</h3>
 *
 * Test system: Pentium 4, Sun Java 1.5_06 -server -Xbatch -Xmx64M <br/>BitSet
 * size = 1,000,000 <br/>Results are java.util.BitSet time divided by OpenBitSet
 * time. <table border="1"> <tr> <th></th> <th>cardinality</th>
 * <th>intersect_count</th> <th>union</th> <th>nextSetBit</th> <th>get</th>
 * <th>iterator</th> </tr> <tr> <th>50% full</th> <td>3.36</td> <td>3.96</td>
 * <td>1.44</td> <td>1.46</td> <td>1.99</td> <td>1.58</td> </tr> <tr> <th>1%
 * full</th> <td>3.31</td> <td>3.90</td> <td>&nbsp;</td> <td>1.04</td>
 * <td>&nbsp;</td> <td>0.99</td> </tr> </table> <br/> Test system: AMD Opteron,
 * 64 bit linux, Sun Java 1.5_06 -server -Xbatch -Xmx64M <br/>BitSet size =
 * 1,000,000 <br/>Results are java.util.BitSet time divided by OpenBitSet time.
 * <table border="1"> <tr> <th></th> <th>cardinality</th>
 * <th>intersect_count</th> <th>union</th> <th>nextSetBit</th> <th>get</th>
 * <th>iterator</th> </tr> <tr> <th>50% full</th> <td>2.50</td> <td>3.50</td>
 * <td>1.00</td> <td>1.03</td> <td>1.12</td> <td>1.25</td> </tr> <tr> <th>1%
 * full</th> <td>2.51</td> <td>3.49</td> <td>&nbsp;</td> <td>1.00</td>
 * <td>&nbsp;</td> <td>1.02</td> </tr> </table>
 *
 * @version $Id$
 */
public class OpenBitSet implements BitSet, Cloneable, Serializable {

    private static final long serialVersionUID = -5197800047652332969L;
    private long[] myBits;
    private int myNumWords;   // number of words (elements) used in the array

    /**
     * Constructs an OpenBitSet large enough to hold numBits.
     *
     * @param numBits
     */
    public OpenBitSet(long numBits) {
        myBits = new long[BitUtil.bits2words(numBits)];
        myNumWords = myBits.length;
    }

    public OpenBitSet() {
        this(64);
    }

    /**
     * Constructs an OpenBitSet from an existing long[]. <br/> The first 64 bits
     * are in long[0], with bit index 0 at the least significant bit, and bit
     * index 63 at the most significant. Given a bit index, the word containing
     * it is long[index/64], and it is at bit number index%64 within that word.
     * <p> numWords are the number of elements in the array that contain set
     * bits (non-zero longs). numWords should be &lt= bits.length, and any
     * existing words in the array at position &gt= numWords should be zero.
     *
     */
    public OpenBitSet(long[] bits, int numWords) {
        if (numWords > bits.length) {
            throw new IllegalArgumentException("OpenBitSet: init: num of words should be less than or equeal to bits length.");
        }
        myBits = bits;
        myNumWords = numWords;
    }
    
    public OpenBitSet(long[] bits) {
        myBits = bits;
        myNumWords = bits.length;
    }
    
    public OpenBitSet(BitSet cloneOBS) {
        myBits = cloneOBS.getBits().clone();
        myNumWords = cloneOBS.getNumWords();
    }

    /**
     * Returns the current capacity in bits (1 greater than the index of the
     * last bit)
     */
    public long capacity() {
        return myBits.length << 6;
    }

    /**
     * Returns the current capacity of this set. Included for compatibility.
     * This is *not* equal to {@link #cardinality}
     */
    public long size() {
        return capacity();
    }

    /**
     * Returns true if there are no set bits
     */
    public boolean isEmpty() {
        return cardinality() == 0;
    }

    /**
     * Expert: returns the long[] storing the bits
     */
    public long[] getBits() {
        return myBits;
    }
    
    public long getBits(int index) {
        return myBits[index];
    }

    /**
     * Expert: sets a new long[] to use as the bit storage
     */
    public void setBits(long[] bits) {
        myBits = bits;
    }
    
    /**
     * Expert: sets specified word with given bits.
     * 
     * @param wordNum word index
     * @param bits bits
     */
    public void setLong(int wordNum, long bits) {
        myBits[wordNum] = bits;
    }

    /**
     * Expert: gets the number of longs in the array that are in use
     */
    public int getNumWords() {
        return myNumWords;
    }

    /**
     * Expert: sets the number of longs in the array that are in use
     */
    public void setNumWords(int nWords) {
        myNumWords = nWords;
    }

    /**
     * Returns true or false for the specified bit index.
     */
    public boolean get(int index) {
        int i = index >> 6;               // div 64
        // signed shift will keep a negative index and force an
        // array-index-out-of-bounds-exception, removing the need for an explicit check.
        if (i >= myBits.length) {
            return false;
        }

        int bit = index & 0x3f;           // mod 64
        long bitmask = 1L << bit;
        return (myBits[i] & bitmask) != 0;
    }

    /**
     * Returns true or false for the specified bit index. The index should be
     * less than the OpenBitSet size
     */
    public boolean fastGet(int index) {
        int i = index >> 6;               // div 64
        // signed shift will keep a negative index and force an
        // array-index-out-of-bounds-exception, removing the need for an explicit check.
        int bit = index & 0x3f;           // mod 64
        long bitmask = 1L << bit;
        return (myBits[i] & bitmask) != 0;
    }

    /**
     * Returns true or false for the specified bit index
     */
    public boolean get(long index) {
        int i = (int) (index >> 6);             // div 64
        if (i >= myBits.length) {
            return false;
        }
        int bit = (int) index & 0x3f;           // mod 64
        long bitmask = 1L << bit;
        return (myBits[i] & bitmask) != 0;
    }

    /**
     * Returns true or false for the specified bit index. The index should be
     * less than the OpenBitSet size.
     */
    public boolean fastGet(long index) {
        int i = (int) (index >> 6);               // div 64
        int bit = (int) index & 0x3f;           // mod 64
        long bitmask = 1L << bit;
        return (myBits[i] & bitmask) != 0;
    }

    /**
     * returns 1 if the bit is set, 0 if not. The index should be less than the
     * OpenBitSet size
     */
    public int getBit(int index) {
        int i = index >> 6;                // div 64
        int bit = index & 0x3f;            // mod 64
        return ((int) (myBits[i] >>> bit)) & 0x01;
    }

    /**
     * sets a bit, expanding the set size if necessary
     */
    public void set(long index) {
        int wordNum = expandingWordNum(index);
        int bit = (int) index & 0x3f;
        long bitmask = 1L << bit;
        myBits[wordNum] |= bitmask;
    }

    /**
     * Sets the bit at the specified index. The index should be less than the
     * OpenBitSet size.
     */
    public void fastSet(int index) {
        int wordNum = index >> 6;      // div 64
        int bit = index & 0x3f;     // mod 64
        long bitmask = 1L << bit;
        myBits[wordNum] |= bitmask;
    }

    /**
     * Sets the bit at the specified index. The index should be less than the
     * OpenBitSet size.
     */
    public void fastSet(long index) {
        int wordNum = (int) (index >> 6);
        int bit = (int) index & 0x3f;
        long bitmask = 1L << bit;
        myBits[wordNum] |= bitmask;
    }

    /**
     * Sets a range of bits, expanding the set size if necessary
     *
     * @param startIndex lower index
     * @param endIndex one-past the last bit to set
     */
    public void set(long startIndex, long endIndex) {
        if (endIndex <= startIndex) {
            return;
        }

        int startWord = (int) (startIndex >> 6);

        // since endIndex is one past the end, this is index of the last
        // word to be changed.
        int endWord = expandingWordNum(endIndex - 1);

        long startmask = -1L << startIndex;
        long endmask = -1L >>> -endIndex;  // 64-(endIndex&0x3f) is the same as -endIndex due to wrap

        if (startWord == endWord) {
            myBits[startWord] |= (startmask & endmask);
            return;
        }

        myBits[startWord] |= startmask;
        Arrays.fill(myBits, startWord + 1, endWord, -1L);
        myBits[endWord] |= endmask;
    }

    protected int expandingWordNum(long index) {
        int wordNum = (int) (index >> 6);
        if (wordNum >= myNumWords) {
            ensureCapacity(index + 1);
            myNumWords = wordNum + 1;
        }
        return wordNum;
    }

    /**
     * clears a bit. The index should be less than the OpenBitSet size.
     */
    public void fastClear(int index) {
        int wordNum = index >> 6;
        int bit = index & 0x03f;
        long bitmask = 1L << bit;
        myBits[wordNum] &= ~bitmask;
        // hmmm, it takes one more instruction to clear than it does to set... any
        // way to work around this?  If there were only 63 bits per word, we could
        // use a right shift of 10111111...111 in binary to position the 0 in the
        // correct place (using sign extension).
        // Could also use Long.rotateRight() or rotateLeft() *if* they were converted
        // by the JVM into a native instruction.
        // bits[word] &= Long.rotateLeft(0xfffffffe,bit);
    }

    /**
     * clears a bit. The index should be less than the OpenBitSet size.
     */
    public void fastClear(long index) {
        int wordNum = (int) (index >> 6); // div 64
        int bit = (int) index & 0x3f;     // mod 64
        long bitmask = 1L << bit;
        myBits[wordNum] &= ~bitmask;
    }

    /**
     * clears a bit, allowing access beyond the current set size without
     * changing the size.
     */
    public void clear(long index) {
        int wordNum = (int) (index >> 6); // div 64
        if (wordNum >= myNumWords) {
            return;
        }
        int bit = (int) index & 0x3f;     // mod 64
        long bitmask = 1L << bit;
        myBits[wordNum] &= ~bitmask;
    }

    /**
     * Clears a range of bits. Clearing past the end does not change the size of
     * the set.
     *
     * @param startIndex lower index
     * @param endIndex one-past the last bit to clear
     */
    public void clear(int startIndex, int endIndex) {
        if (endIndex <= startIndex) {
            return;
        }

        int startWord = (startIndex >> 6);
        if (startWord >= myNumWords) {
            return;
        }

        // since endIndex is one past the end, this is index of the last
        // word to be changed.
        int endWord = ((endIndex - 1) >> 6);

        long startmask = -1L << startIndex;
        long endmask = -1L >>> -endIndex;  // 64-(endIndex&0x3f) is the same as -endIndex due to wrap

        // invert masks since we are clearing
        startmask = ~startmask;
        endmask = ~endmask;

        if (startWord == endWord) {
            myBits[startWord] &= (startmask | endmask);
            return;
        }

        myBits[startWord] &= startmask;

        int middle = Math.min(myNumWords, endWord);
        Arrays.fill(myBits, startWord + 1, middle, 0L);
        if (endWord < myNumWords) {
            myBits[endWord] &= endmask;
        }
    }

    /**
     * Clears a range of bits. Clearing past the end does not change the size of
     * the set.
     *
     * @param startIndex lower index
     * @param endIndex one-past the last bit to clear
     */
    public void clear(long startIndex, long endIndex) {
        if (endIndex <= startIndex) {
            return;
        }

        int startWord = (int) (startIndex >> 6);
        if (startWord >= myNumWords) {
            return;
        }

        // since endIndex is one past the end, this is index of the last
        // word to be changed.
        int endWord = (int) ((endIndex - 1) >> 6);

        long startmask = -1L << startIndex;
        long endmask = -1L >>> -endIndex;  // 64-(endIndex&0x3f) is the same as -endIndex due to wrap

        // invert masks since we are clearing
        startmask = ~startmask;
        endmask = ~endmask;

        if (startWord == endWord) {
            myBits[startWord] &= (startmask | endmask);
            return;
        }

        myBits[startWord] &= startmask;

        int middle = Math.min(myNumWords, endWord);
        Arrays.fill(myBits, startWord + 1, middle, 0L);
        if (endWord < myNumWords) {
            myBits[endWord] &= endmask;
        }
    }

    /**
     * Sets a bit and returns the previous value. The index should be less than
     * the OpenBitSet size.
     */
    public boolean getAndSet(int index) {
        int wordNum = index >> 6;      // div 64
        int bit = index & 0x3f;     // mod 64
        long bitmask = 1L << bit;
        boolean val = (myBits[wordNum] & bitmask) != 0;
        myBits[wordNum] |= bitmask;
        return val;
    }

    /**
     * Sets a bit and returns the previous value. The index should be less than
     * the OpenBitSet size.
     */
    public boolean getAndSet(long index) {
        int wordNum = (int) (index >> 6);      // div 64
        int bit = (int) index & 0x3f;     // mod 64
        long bitmask = 1L << bit;
        boolean val = (myBits[wordNum] & bitmask) != 0;
        myBits[wordNum] |= bitmask;
        return val;
    }

    /**
     * flips a bit. The index should be less than the OpenBitSet size.
     */
    public void fastFlip(int index) {
        int wordNum = index >> 6;      // div 64
        int bit = index & 0x3f;     // mod 64
        long bitmask = 1L << bit;
        myBits[wordNum] ^= bitmask;
    }

    /**
     * flips a bit. The index should be less than the OpenBitSet size.
     */
    public void fastFlip(long index) {
        int wordNum = (int) (index >> 6);   // div 64
        int bit = (int) index & 0x3f;       // mod 64
        long bitmask = 1L << bit;
        myBits[wordNum] ^= bitmask;
    }

    /**
     * flips a bit, expanding the set size if necessary
     */
    public void flip(long index) {
        int wordNum = expandingWordNum(index);
        int bit = (int) index & 0x3f;       // mod 64
        long bitmask = 1L << bit;
        myBits[wordNum] ^= bitmask;
    }

    /**
     * flips a bit and returns the resulting bit value. The index should be less
     * than the OpenBitSet size.
     */
    public boolean flipAndGet(int index) {
        int wordNum = index >> 6;      // div 64
        int bit = index & 0x3f;     // mod 64
        long bitmask = 1L << bit;
        myBits[wordNum] ^= bitmask;
        return (myBits[wordNum] & bitmask) != 0;
    }

    /**
     * flips a bit and returns the resulting bit value. The index should be less
     * than the OpenBitSet size.
     */
    public boolean flipAndGet(long index) {
        int wordNum = (int) (index >> 6);   // div 64
        int bit = (int) index & 0x3f;       // mod 64
        long bitmask = 1L << bit;
        myBits[wordNum] ^= bitmask;
        return (myBits[wordNum] & bitmask) != 0;
    }

    /**
     * Flips a range of bits, expanding the set size if necessary
     *
     * @param startIndex lower index
     * @param endIndex one-past the last bit to flip
     */
    public void flip(long startIndex, long endIndex) {
        if (endIndex <= startIndex) {
            return;
        }
        int startWord = (int) (startIndex >> 6);

        // since endIndex is one past the end, this is index of the last
        // word to be changed.
        int endWord = expandingWordNum(endIndex - 1);

        /**
         * * Grrr, java shifting wraps around so -1L>>>64 == -1 for that
         * reason, make sure not to use endmask if the bits to flip will be zero
         * in the last word (redefine endWord to be the last changed...) long
         * startmask = -1L << (startIndex & 0x3f); // example: 11111...111000
         * long endmask = -1L >>> (64-(endIndex & 0x3f)); // example:
         * 00111...111111 *
         */
        long startmask = -1L << startIndex;
        long endmask = -1L >>> -endIndex;  // 64-(endIndex&0x3f) is the same as -endIndex due to wrap

        if (startWord == endWord) {
            myBits[startWord] ^= (startmask & endmask);
            return;
        }

        myBits[startWord] ^= startmask;

        for (int i = startWord + 1; i < endWord; i++) {
            myBits[i] = ~myBits[i];
        }

        myBits[endWord] ^= endmask;
    }

    /**
     * @return the number of set bits
     */
    public long cardinality() {
        return BitUtil.pop_array(myBits, 0, myNumWords);
    }
    
    public long cardinality(int index) {
        return BitUtil.pop_array_to_index(myBits, index);
    }

    /**
     * Returns the popcount or cardinality of the intersection of the two sets.
     * Neither set is modified.
     */
    public static long intersectionCount(BitSet a, BitSet b) {
        return BitUtil.pop_intersect(a.getBits(), b.getBits(), 0, Math.min(a.getNumWords(), b.getNumWords()));
    }

    /**
     * Returns the popcount or cardinality of the union of the two sets. Neither
     * set is modified.
     */
    public static long unionCount(BitSet a, BitSet b) {
        long tot = BitUtil.pop_union(a.getBits(), b.getBits(), 0, Math.min(a.getNumWords(), b.getNumWords()));
        if (a.getNumWords() < b.getNumWords()) {
            tot += BitUtil.pop_array(b.getBits(), a.getNumWords(), b.getNumWords() - a.getNumWords());
        } else if (a.getNumWords() > b.getNumWords()) {
            tot += BitUtil.pop_array(a.getBits(), b.getNumWords(), a.getNumWords() - b.getNumWords());
        }
        return tot;
    }

    /**
     * Returns the popcount or cardinality of "a and not b" or "intersection(a,
     * not(b))". Neither set is modified.
     */
    public static long andNotCount(BitSet a, BitSet b) {
        long tot = BitUtil.pop_andnot(a.getBits(), b.getBits(), 0, Math.min(a.getNumWords(), b.getNumWords()));
        if (a.getNumWords() > b.getNumWords()) {
            tot += BitUtil.pop_array(a.getBits(), b.getNumWords(), a.getNumWords() - b.getNumWords());
        }
        return tot;
    }

    /**
     * Returns the popcount or cardinality of the exclusive-or of the two sets.
     * Neither set is modified.
     */
    public static long xorCount(BitSet a, BitSet b) {
        long tot = BitUtil.pop_xor(a.getBits(), b.getBits(), 0, Math.min(a.getNumWords(), b.getNumWords()));
        if (a.getNumWords() < b.getNumWords()) {
            tot += BitUtil.pop_array(b.getBits(), a.getNumWords(), b.getNumWords() - a.getNumWords());
        } else if (a.getNumWords() > b.getNumWords()) {
            tot += BitUtil.pop_array(a.getBits(), b.getNumWords(), a.getNumWords() - b.getNumWords());
        }
        return tot;
    }

    /**
     * Returns the index of the first set bit starting at the index specified.
     * -1 is returned if there are no more set bits.
     */
    public int nextSetBit(int index) {
        int i = index >> 6;
        if (i >= myNumWords) {
            return -1;
        }
        int subIndex = index & 0x3f;      // index within the word
        long word = myBits[i] >> subIndex;  // skip all the bits to the right of index

        if (word != 0) {
            return (i << 6) + subIndex + BitUtil.ntz(word);
        }

        while (++i < myNumWords) {
            word = myBits[i];
            if (word != 0) {
                return (i << 6) + BitUtil.ntz(word);
            }
        }

        return -1;
    }

    /**
     * Returns the index of the first set bit starting at the index specified.
     * -1 is returned if there are no more set bits.
     */
    public long nextSetBit(long index) {
        int i = (int) (index >>> 6);
        if (i >= myNumWords) {
            return -1;
        }
        int subIndex = (int) index & 0x3f; // index within the word
        long word = myBits[i] >>> subIndex;  // skip all the bits to the right of index

        if (word != 0) {
            return (((long) i) << 6) + (subIndex + BitUtil.ntz(word));
        }

        while (++i < myNumWords) {
            word = myBits[i];
            if (word != 0) {
                return (((long) i) << 6) + BitUtil.ntz(word);
            }
        }

        return -1;
    }

	@Override
	public int previousSetBit(int index) {
        int i = index >> 6;
        if (i >= myNumWords) {
            return -1;
        }
        int subIndex = index & 0x3f;      // index within the word
        long word = myBits[i] << (63 - subIndex);  // skip all the bits to the left of index
        if (word != 0) {
        	int prevIndex = subIndex;
        	while(word > 0) {
        		word = word << 1;
        		prevIndex--;
        	}
            return (i << 6) + prevIndex;
        }

        while (--i >= 0) {
        	word = myBits[i];
            if (word != 0) {
            	int prevIndex = 63;
            	while(word > 0) {
            		word = word << 1;
            		prevIndex--;
            	}
                return (i << 6) + prevIndex;
            }
        }

        return -1;
	}

	@Override
	public long previousSetBit(long index) {
        int i = (int) (index >> 6);
        if (i >= myNumWords) {
            return -1;
        }
        int subIndex = (int) (index & 0x3f);      // index within the word
        long word = myBits[i] << (63 - subIndex);  // skip all the bits to the left of index
        if (word != 0) {
        	int prevIndex = subIndex;
        	while(word > 0) {
        		word = word << 1;
        		prevIndex--;
        	}
            return (i << 6) + prevIndex;
        }

        while (--i >= 0) {
        	word = myBits[i];
            if (word != 0) {
            	int prevIndex = 63;
            	while(word > 0) {
            		word = word << 1;
            		prevIndex--;
            	}
                return (i << 6) + prevIndex;
            }
        }

        return -1;
      }

	@Override
    public Object clone() {
        long[] bits = (long[]) getBits().clone();
        return new OpenBitSet(bits, getNumWords());
    }

    /**
     * this = this AND other
     */
    public void intersect(BitSet other) {
        int newLen = Math.min(myNumWords, other.getNumWords());
        long[] thisArr = myBits;
        long[] otherArr = other.getBits();
        // testing against zero can be more efficient
        int pos = newLen;
        while (--pos >= 0) {
            thisArr[pos] &= otherArr[pos];
        }
        if (myNumWords > newLen) {
            // fill zeros from the new shorter length to the old length
            Arrays.fill(myBits, newLen, myNumWords, 0);
        }
        myNumWords = newLen;
    }

    /**
     * this = this OR other
     */
    public void union(BitSet other) {
        int newLen = Math.max(myNumWords, other.getNumWords());
        ensureCapacityWords(newLen);

        long[] thisArr = myBits;
        long[] otherArr = other.getBits();
        int pos = Math.min(myNumWords, other.getNumWords());
        while (--pos >= 0) {
            thisArr[pos] |= otherArr[pos];
        }
        if (myNumWords < newLen) {
            System.arraycopy(otherArr, myNumWords, thisArr, myNumWords, newLen - myNumWords);
        }
        myNumWords = newLen;
    }

    /**
     * Remove all elements set in other. this = this AND_NOT other
     */
    public void remove(BitSet other) {
        int idx = Math.min(myNumWords, other.getNumWords());
        long[] thisArr = myBits;
        long[] otherArr = other.getBits();
        while (--idx >= 0) {
            thisArr[idx] &= ~otherArr[idx];
        }
    }

    /**
     * this = this XOR other
     */
    public void xor(BitSet other) {
        int newLen = Math.max(myNumWords, other.getNumWords());
        ensureCapacityWords(newLen);

        long[] thisArr = myBits;
        long[] otherArr = other.getBits();
        int pos = Math.min(myNumWords, other.getNumWords());
        while (--pos >= 0) {
            thisArr[pos] ^= otherArr[pos];
        }
        if (myNumWords < newLen) {
            System.arraycopy(otherArr, myNumWords, thisArr, myNumWords, newLen - myNumWords);
        }
        myNumWords = newLen;
    }

    /**
     * this = not(this XOR other)
     */
    public void notXor(BitSet other) {
        int newLen = Math.max(myNumWords, other.getNumWords());
        ensureCapacityWords(newLen);

        long[] thisArr = myBits;
        long[] otherArr = other.getBits();
        int pos = Math.min(myNumWords, other.getNumWords());
        while (--pos >= 0) {
            thisArr[pos] = ~(thisArr[pos] ^ otherArr[pos]);
        }
        if (myNumWords < newLen) {
            System.arraycopy(otherArr, myNumWords, thisArr, myNumWords, newLen - myNumWords);
        }
        myNumWords = newLen;
    }

    /**
     * this = not(this)
     */
    public void not() {
        long[] thisArr = myBits;
        int pos = myNumWords;
        while (--pos >= 0) {
            thisArr[pos] = ~(thisArr[pos]);
        }
    }

    // some BitSet compatability methods
    /**
     * see {@link intersect}
     */
    public void and(BitSet other) {
        intersect(other);
    }

    /**
     * see {@link union}
     */
    public void or(BitSet other) {
        union(other);
    }

    /**
     * see {@link andNot}
     */
    public void andNot(BitSet other) {
        remove(other);
    }

    /**
     * returns true if the sets have any elements in common
     */
    public boolean intersects(BitSet other) {
        int pos = Math.min(myNumWords, other.getNumWords());
        long[] thisArr = myBits;
        long[] otherArr = other.getBits();
        while (--pos >= 0) {
            if ((thisArr[pos] & otherArr[pos]) != 0) {
                return true;
            }
        }
        return false;
    }

    /**
     * Expand the long[] with the size given as a number of words (64 bit
     * longs). getNumWords() is unchanged by this call.
     */
    public void ensureCapacityWords(int numWords) {
        if (myBits.length < numWords) {
            long[] newBits = new long[numWords];
            System.arraycopy(myBits, 0, newBits, 0, myNumWords);
            myBits = newBits;
        }
    }

    /**
     * Ensure that the long[] is big enough to hold numBits, expanding it if
     * necessary. getNumWords() is unchanged by this call.
     */
    public void ensureCapacity(long numBits) {
        ensureCapacityWords(BitUtil.bits2words(numBits));
    }

    /**
     * Lowers numWords, the number of words in use, by checking for trailing
     * zero words.
     */
    public void trimTrailingZeros() {
        int idx = myNumWords - 1;
        while (idx >= 0 && myBits[idx] == 0) {
            idx--;
        }
        myNumWords = idx + 1;
    }

    public int indexOfNthSetBit(int n) {
        
        if (n < 1) {
            return -1;
        }

        int result = 0;
        int count = 0;
        for (int i = 0; i < myNumWords; i++) {
            int currentCount = BitUtil.pop(myBits[i]);
            if ((count + currentCount) >= n) {
                long bitmask = 1L;
                while (true) {
                    if ((myBits[i] & bitmask) != 0) {
                        count++;
                        if (count == n) {
                            return result;
                        }
                    }
                    bitmask = bitmask << 1;
                    result++;
                }
            } else {
                count = count + currentCount;
                result = result + 64;
            }
        }

        return -1;

    }

    public int[] getIndicesOfSetBits() {

        int[] result = new int[(int) cardinality()];
        int count = 0;

        for (int i = 0; i < myNumWords; i++) {
            long bitmask = 1L;
            for (int j = 0; j < 64; j++) {
                if ((myBits[i] & bitmask) != 0) {
                    result[count++] = i * 64 + j;
                }
                bitmask = bitmask << 1;
            }
        }

        return result;

    }

    /**
     * returns true if both sets have the same bits set
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof BitSet)) {
            return false;
        }
        BitSet a;
        BitSet b = (BitSet) o;
        // make a the larger set.
        if (b.getNumWords() > myNumWords) {
            a = b;
            b = this;
        } else {
            a = this;
        }

        // check for any set bits out of the range of b
        for (int i = a.getNumWords() - 1; i >= b.getNumWords(); i--) {
            if (a.getBits()[i] != 0) {
                return false;
            }
        }

        for (int i = b.getNumWords() - 1; i >= 0; i--) {
            if (a.getBits()[i] != b.getBits()[i]) {
                return false;
            }
        }

        return true;
    }

    @Override
    public int hashCode() {
        // Start with a zero hash and use a mix that results in zero if the input is zero.
        // This effectively truncates trailing zeros without an explicit check.
        long h = 0;
        for (int i = myBits.length; --i >= 0;) {
            h ^= myBits[i];
            h = (h << 1) | (h >>> 63); // rotate left
        }
        // fold leftmost bits into right and add a constant to prevent
        // empty sets from returning 0, which is too common.
        return (int) ((h >> 32) ^ h) + 0x98761234;
    }

    @Override
    public long[] getBits(int startWord, int endWord) {
        int length=endWord-startWord+1;
        long[] sL = new long[length];
        System.arraycopy(myBits, startWord, sL, 0, length);
        return sL;
    }
}
