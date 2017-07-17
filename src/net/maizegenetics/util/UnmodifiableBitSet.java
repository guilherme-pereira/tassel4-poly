/*
 * UnmodifiableBitSet
 */
package net.maizegenetics.util;

/**
 *
 * @author terry
 */
public class UnmodifiableBitSet implements BitSet {
    
    private final BitSet myBitSet;
    
    private UnmodifiableBitSet(BitSet bitSet) {
        myBitSet = bitSet;
    }
    
    public static UnmodifiableBitSet getInstance(BitSet bitSet) {
        if (bitSet instanceof UnmodifiableBitSet) {
            return (UnmodifiableBitSet) bitSet;
        } else {
            return new UnmodifiableBitSet(bitSet);
        }
    }
    
    @Override
    public long capacity() {
        return myBitSet.capacity();
    }
    
    @Override
    public long size() {
        return myBitSet.size();
    }
    
    @Override
    public boolean isEmpty() {
        return myBitSet.isEmpty();
    }
    
    @Override
    public long[] getBits() {
        return myBitSet.getBits().clone();
    }
    
    @Override
    public long[] getBits(int startWord, int endWord) {
        int length=endWord-startWord+1;
        long[] sL = new long[length];
        System.arraycopy(myBitSet.getBits(), startWord, sL, 0, length);
        return sL;
    }
    
    @Override
    public long getBits(int index) {
        return myBitSet.getBits(index);
    }
    
    @Override
    public void setBits(long[] bits) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void setLong(int wordNum, long bits) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public int getNumWords() {
        return myBitSet.getNumWords();
    }
    
    @Override
    public void setNumWords(int nWords) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public boolean get(int index) {
        return myBitSet.get(index);
    }
    
    @Override
    public boolean fastGet(int index) {
        return myBitSet.fastGet(index);
    }
    
    @Override
    public boolean get(long index) {
        return myBitSet.get(index);
    }
    
    @Override
    public boolean fastGet(long index) {
        return myBitSet.fastGet(index);
    }
    
    @Override
    public int getBit(int index) {
        return myBitSet.getBit(index);
    }
    
    @Override
    public void set(long index) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void fastSet(int index) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void fastSet(long index) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void set(long startIndex, long endIndex) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void fastClear(int index) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void fastClear(long index) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void clear(long index) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void clear(int startIndex, int endIndex) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void clear(long startIndex, long endIndex) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public boolean getAndSet(int index) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public boolean getAndSet(long index) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void fastFlip(int index) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void fastFlip(long index) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void flip(long index) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public boolean flipAndGet(int index) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public boolean flipAndGet(long index) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void flip(long startIndex, long endIndex) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public long cardinality() {
        return myBitSet.cardinality();
    }
    
    @Override
    public long cardinality(int index) {
        return myBitSet.cardinality(index);
    }
    
    @Override
    public int nextSetBit(int index) {
        return myBitSet.nextSetBit(index);
    }
    
    @Override
    public long nextSetBit(long index) {
        return myBitSet.nextSetBit(index);
    }
    
    @Override
	public int previousSetBit(int index) {
		return myBitSet.previousSetBit(index);
	}

	@Override
	public long previousSetBit(long index) {
		return myBitSet.previousSetBit(index);
	}

	@Override
    public void intersect(BitSet other) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void union(BitSet other) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void remove(BitSet other) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void xor(BitSet other) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void and(BitSet other) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void or(BitSet other) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void andNot(BitSet other) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public boolean intersects(BitSet other) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void ensureCapacityWords(int numWords) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void ensureCapacity(long numBits) {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public void trimTrailingZeros() {
        throw new UnsupportedOperationException("UnmodifiableBitSet.");
    }
    
    @Override
    public int indexOfNthSetBit(int n) {
        return myBitSet.indexOfNthSetBit(n);
    }
    
    @Override
    public int[] getIndicesOfSetBits() {
        return myBitSet.getIndicesOfSetBits();
    }
    
    @Override
    public boolean equals(Object o) {
        return myBitSet.equals(o);
    }
}
