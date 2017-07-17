/*
 * PolymorphismFinder
 */
package net.maizegenetics.gbs.pipeline;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;

import java.util.ArrayList;
import java.util.Arrays;

import net.maizegenetics.gbs.util.ReadsByTaxa;

/**
 *
 * @author fl262
 */
public class PolymorphismFinder {

    ReadsByTaxa rbt;
    long[] lookup;
    int[] indexOfRbt;

    public PolymorphismFinder(ReadsByTaxa rbt) {
        this.rbt = rbt;
        initialize();
    }

    public void initialize() {
        lookup = new long[rbt.haplotypeNum * 2];
        indexOfRbt = new int[rbt.haplotypeNum * 2];
        for (int i = 0; i < rbt.haplotypeNum; i++) {
            for (int j = 0; j < 2; j++) {
                int dex = 2 * i + j;
                lookup[dex] = rbt.haplotype[j][i];
                indexOfRbt[dex] = i;
            }
        }
        GenericSorting.quickSort(0, lookup.length, comp, swapper);
        this.reduceDuplicates();
        System.out.println("Initialization is done");
    }

    public ArrayList<Integer> findOneMismatch(long[] queryLongSeq) {
        ArrayList<Integer> hitIndex = new ArrayList();
        for (int i = 0; i < queryLongSeq.length; i++) {
            int hit = Arrays.binarySearch(lookup, queryLongSeq[i]);
            if (hit < 0) {
                continue;
            }
            while ((hit > 0) && (queryLongSeq[i] == lookup[hit - 1])) {
                hit--;
            }
            while ((hit < lookup.length) && (lookup[hit] == queryLongSeq[i])) {
                int count = getMismatchInLong(queryLongSeq[0], rbt.haplotype[0][indexOfRbt[hit]]) + getMismatchInLong(queryLongSeq[1], rbt.haplotype[1][indexOfRbt[hit]]);
                if (count == 1) {
                    hitIndex.add(indexOfRbt[hit]);
                }
                hit++;
            }
        }
        return hitIndex;
    }

    public byte getMismatchInLong(long longSeq1, long longSeq2) {
        long mask = 3;
        long diff = longSeq1 ^ longSeq2;
        byte count = 0;
        for (int i = 0; i < 32; i++) {
            if ((diff & mask) > 0) {
                count++;
            }
            diff = diff >> 2;
        }
        return count;
    }

    private void reduceDuplicates() {
        int start = 0, end = -1, duplicated = 0;
        long currHap = lookup[0];
        for (int i = 0; i < lookup.length; i++) {
            if (lookup[i] == currHap) {
                end = i;
            } else {
                if (((end - start) > 1000)) {
                    //System.out.println(BaseEncoder.getSequenceFromInt(currHap)+" "+(end-start));
                    for (int j = start; j <= end; j++) {
                        lookup[j] = Long.MAX_VALUE;
                        duplicated++;
                    }
                }
                currHap = lookup[i];
                start = end = i;
            }
        }
        GenericSorting.quickSort(0, lookup.length, comp, swapper);
        long[] newlookup = new long[lookup.length - duplicated];
        int[] newindexOfRbt = new int[lookup.length - duplicated];
        System.arraycopy(lookup, 0, newlookup, 0, lookup.length - duplicated);
        System.arraycopy(indexOfRbt, 0, newindexOfRbt, 0, indexOfRbt.length - duplicated);
        System.out.println("Old Lookup Size:" + lookup.length + "  new size:" + newlookup.length);
        lookup = newlookup;
        indexOfRbt = newindexOfRbt;
    }
    Swapper swapper = new Swapper() {

        public void swap(int a, int b) {
            long tl;
            int ti, tb;
            tl = lookup[a];
            lookup[a] = lookup[b];
            lookup[b] = tl;
            ti = indexOfRbt[a];
            indexOfRbt[a] = indexOfRbt[b];
            indexOfRbt[b] = ti;
        }
    };
    IntComparator comp = new IntComparator() {

        public int compare(int a, int b) {
            if (lookup[a] < lookup[b]) {
                return -1;
            } else if (lookup[a] > lookup[b]) {
                return 1;
            } else {
                return 0;
            }
        }
    };
}
