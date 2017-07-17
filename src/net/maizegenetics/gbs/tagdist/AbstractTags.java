/*
 * AbstractTags
 */
package net.maizegenetics.gbs.tagdist;

import cern.colt.GenericSorting;
import java.util.Arrays;

/**
 * Basic methods for working with Tags, including sorting and search.
 * @author edbuckler
 */
public abstract class AbstractTags implements Tags {

    protected int tagLengthInLong;  //TODO fully implement on reading
    protected long[][] tags;  // for memory efficiency the rows first and second half of the read
    // columns are the index of the reads.
    protected byte[] tagLength;  // length of tag (number of bases)  // 1 byte

    @Override
    public boolean areTagsUnique() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public long[] getTag(int index) {
        long[] theTag = new long[tagLengthInLong];
        for (int i = 0; i < tagLengthInLong; i++) {
            theTag[i] = tags[i][index];
        }
        return theTag;
    }

    @Override
    public int getTagCount() {
        return tags[0].length;
    }

    @Override
    public int getTagIndex(long[] read) {
        //code inspired by COLT lower bound function
        int first = 0;
        int len = tags[0].length - first;
        int comp = 0;
        while (len > 0) {
            int half = len / 2;
            int middle = first + half;
            if ((comp = compareTags(middle, read)) < 0) {
                first = middle + 1;
                len -= half + 1;
            } else {
                len = half;
            }
        }
        if ((first < tags[0].length) && (compareTags(first, read) == 0)) {
            return first;
        }
        return -(first + 1);
    }

    public int getTagIndexFirst(long readFirst) {
        return Arrays.binarySearch(tags[0], readFirst);
    }

    public boolean areTagsEqual(int index1, int index2) {
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tags[i][index1] != tags[i][index2]) {
                return false;
            }
        }
        return true;
    }

    public int compareTags(int index1, int index2) {
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tags[i][index1] < tags[i][index2]) {
                return -1;
            }
            if (tags[i][index1] > tags[i][index2]) {
                return 1;
            }
        }
        return 0;
    }

    public int compareTags(int index1, long[] t2) {
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tags[i][index1] < t2[i]) {
                return -1;
            }
            if (tags[i][index1] > t2[i]) {
                return 1;
            }
        }
        return 0;
    }

    public static int compareTags(long[] t1, long[] t2) {
        for (int i = 0; i < t1.length; i++) {
            if (t1[i] < t2[i]) {
                return -1;
            }
            if (t1[i] > t2[i]) {
                return 1;
            }
        }
        return 0;
    }

    @Override
    public int[] getTagIndexSet(long[] read) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getTagLength(int index) {
        return tagLength[index];
    }

    public byte[] getTagLength() {
        return tagLength;
    }

    @Override
    public int getTagSizeInLong() {
        return tagLengthInLong;
    }

    @Override
    public String getNullTag() {
        char[] nullBases = new char[32 * tagLengthInLong];
        Arrays.fill(nullBases, 'A');
        return new String(nullBases);
    }

    @Override
    public void swap(int index1, int index2) {
        long temp;
        for (int i = 0; i < tagLengthInLong; i++) {
            temp = tags[i][index1];
            tags[i][index1] = tags[i][index2];
            tags[i][index2] = temp;
        }
        byte tl;
        tl = tagLength[index1];
        tagLength[index1] = tagLength[index2];
        tagLength[index2] = tl;
    }

    @Override
    public int compare(int index1, int index2) {
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tags[i][index1] < tags[i][index2]) {
                return -1;
            }
            if (tags[i][index1] > tags[i][index2]) {
                return 1;
            }
        }
        return 0;
    }

    public void sort() {
        System.out.println("Position index sort begin.");
        GenericSorting.quickSort(0, this.getTagCount(), this, this);
        System.out.println("Position index sort end.");
    }
    /**
     * Generically searches the list for the specified value using
     * the binary search algorithm.  The list must <strong>must</strong> be
     * sorted (as by the sort method) prior to making this call.  If
     * it is not sorted, the results are undefined: in particular, the call
     * may enter an infinite loop.  If the list contains multiple elements
     * equal to the specified key, there is no guarantee which of the multiple elements
     * will be found.
     *
     * @param list the list to be searched.
     * @param key the value to be searched for.
     * @param from the leftmost search position, inclusive.
     * @param to the rightmost search position, inclusive.
     * @return index of the search key, if it is contained in the list;
     *	       otherwise, <tt>(-(<i>insertion point</i>) - 1)</tt>.  The <i>insertion
     *	       point</i> is defined as the the point at which the value would
     * 	       be inserted into the list: the index of the first
     *	       element greater than the key, or <tt>list.length</tt>, if all
     *	       elements in the list are less than the specified key.  Note
     *	       that this guarantees that the return value will be &gt;= 0 if
     *	       and only if the key is found.
     * @see java.util.Arrays
     */
//    public static int binarySearchFromTo(int from, int to, IntComparator comp) {
//            final int dummy = 0;
//            while (from <= to) {
//                    int mid = (from + to) / 2;
//                    int comparison = comp.compare(dummy,mid);
//                    if (comparison < 0) from = mid + 1;
//                    else if (comparison > 0) to = mid - 1;
//                    else return mid; // key found
//            }
//            return -(from + 1);  // key not found.
//    }
//
//    private static int lower_bound(int[] array, int first, int last, int x) {
//                    int len = last - first;
//                    while (len > 0) {
//                            int half = len / 2;
//                            int middle = first + half;
//                            if (array[middle] < x) {
//                                    first = middle + 1;
//                                    len -= half + 1;
//                            } else
//                                    len = half;
//                    }
//                    return first;
//            }
}
