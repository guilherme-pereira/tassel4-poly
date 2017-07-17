package net.maizegenetics.gbs.homology;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;

import java.util.Arrays;
import java.util.Map.Entry;
import java.util.TreeMap;

import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.gbs.util.BaseEncoder;

/**
 * Very simple but fast homology finder.  Very similar to BLAT, with long word searches.
 * A in memory index is created from any set of Tags, and this list many then be
 * quickly queried for high homology hits.
 *<p>
 * With the current code no indels will be discovered, however, this could be easily changed.
 * Currently, tuned to a word length of 16.  If you are going to use this please 
 * read this code.
 * <p>
 * The lookup match tables are based on 2-bit encoded long sequences, and the query also needs to
 * in a 2-bit long
 * 
 * @author Ed Buckler
 */
public class TagMatchFinder {

    Tags refTags;
    private static int[][] lookups;
    public static final int wordLength = 16;        //best results 8 for sensitivity
    public static final int maxDivergence = 15;
    private SmithWaterman pa;
//    private Thread liveThread = null; //use to keep this alive while another thread is alive

    public int[][] getTagLookTable() {
        return lookups;
    }

    public TagMatchFinder(Tags theTags) {
        this.refTags = theTags;
        init();
        pa = new SmithWaterman(64, 32);
        long[] query = refTags.getTag(10);
        TreeMap<Integer, Integer> al = findMatchesWithIntLengthWords(query, 10, false);
        System.out.println(al.toString());
    }

    public static void main(String[] args) {
//       TagCounts refTags=new TagCounts("C:/EdStuff/Solexa/NGG_IBM/counts/countMerge091223_20.txt",true);
//       TagMatchFinder tmf=new TagMatchFinder(refTags);
    }

    private void init() {
        lookups = new int[2][2 * refTags.getTagCount() * refTags.getTagSizeInLong()];  //4 x 16base words
        int count = 0;
        System.out.println("Beginning to create index");
        for (int i = 0; i < refTags.getTagCount(); i++) {
            long[] cTag = refTags.getTag(i);
            for (int j = 0; j < refTags.getTagSizeInLong(); j++) {
                int word[] = BaseEncoder.getIntFromLong(cTag[j]);
                lookups[0][count] = word[0];
                lookups[1][count] = i;
                count++;
                lookups[0][count] = word[1];
                lookups[1][count] = i;
                count++;
            }
        }
        System.out.println("Sorting index");
        GenericSorting.quickSort(0, lookups[0].length, comp, swapper);
        System.out.println("Beginning to reduce duplicates");
        reduceDuplicates();
        System.out.println("Duplicates removed");
        pa = new SmithWaterman(32, 32);
    }

    private void reduceDuplicates() {
        int start = 0, end = -1, currHap = lookups[0][0], duplicated = 0;
        System.out.println(BaseEncoder.getSequenceFromInt(currHap) + " " + currHap);
        for (int i = 0; i < lookups[0].length; i++) {
            if (lookups[0][i] == currHap) {
                end = i;
            } else {
                if (((end - start) > 1000)) {
                    //System.out.println(BaseEncoder.getSequenceFromInt(currHap)+" "+(end-start));
                    for (int j = start; j <= end; j++) {
                        lookups[0][j] = Integer.MAX_VALUE;
                        duplicated++;
                    }
                }
                currHap = lookups[0][i];
                start = end = i;
            }
        }
        GenericSorting.quickSort(0, lookups[0].length, comp, swapper);
        int[][] newlookups = new int[2][lookups[0].length - duplicated];
        System.arraycopy(lookups[0], 0, newlookups[0], 0, lookups[0].length - duplicated);
        System.arraycopy(lookups[1], 0, newlookups[1], 0, lookups[1].length - duplicated);
        System.out.println("Old Lookup Size:" + lookups[0].length + "  new size:" + newlookups[0].length);
        lookups = newlookups;
    }

    /**
     * Return a TreeMap good hits based on a sequence query.  The returned tree map
     * is the list of tag indices as key, divergence as value.  
     * Ed- It seems like the key & value should perhaps be reversed.
     * @param query array of 2-bit encoded long query
     * @param maxDiv maximum divergence to look for
     * @param keepOnlyBest result only includes the best result
     * @return hit map (tagIndex, divergence)
     */
    public TreeMap<Integer, Integer> findMatchesWithIntLengthWords(long[] query, int maxDiv, boolean keepOnlyBest) {
        TreeMap<Integer, Integer> hitsAndDiv = new TreeMap<Integer, Integer>();
        int bestDiv = maxDiv;
        //   OpenBitSet bQuery=new OpenBitSet(query);
        for (int i = 0; i < query.length; i++) {
            int words[] = BaseEncoder.getIntFromLong(query[i]);
            for (int word : words) {
                //      System.out.println("word"+word+" string:"+BaseEncoder.getSequenceFromInt(word));
                int hit = Arrays.binarySearch(lookups[0], word);
                if (hit < 0) {
                    continue;  //no hit for this word
                }
                while ((hit > 0) && (word == lookups[0][hit - 1])) {
                    hit--;
                }  //backup to the first hit
                while ((hit < lookups[0].length) && (lookups[0][hit] == word)) {
                    int indexPossHapHit = lookups[1][hit];
                    if (!hitsAndDiv.containsKey(indexPossHapHit)) {
                        int div = BaseEncoder.seqDifferences(refTags.getTag(indexPossHapHit)[0], query[0], maxDiv)
                                + BaseEncoder.seqDifferences(refTags.getTag(indexPossHapHit)[1], query[1], maxDiv);
                        if (div <= maxDiv) {
                            hitsAndDiv.put(indexPossHapHit, div); //System.out.println("Hit"+indexPossHapHit+" div:"+div);
                        }
                        if (div < bestDiv) {
                            bestDiv = div;
                        }
                    }
                    hit++;
                }
            }
        }
        if (keepOnlyBest && hitsAndDiv.size() > 1) {
            TreeMap<Integer, Integer> bestHitsAndDiv = new TreeMap<Integer, Integer>();
            for (Entry<Integer, Integer> e : hitsAndDiv.entrySet()) {
                if (e.getValue() == bestDiv) {
                    bestHitsAndDiv.put(e.getKey(), e.getValue());
                }
            }
            return bestHitsAndDiv;
        } else {
            return hitsAndDiv;
        }
        //1bp indels can be dealt with by shift
    }
    Swapper swapper = new Swapper() {

        public void swap(int a, int b) {
            int t1, t2;
            t1 = lookups[0][a];
            lookups[0][a] = lookups[0][b];
            lookups[0][b] = t1;
            t2 = lookups[1][a];
            lookups[1][a] = lookups[1][b];
            lookups[1][b] = t2;
        }
    };
    IntComparator comp = new IntComparator() {

        public int compare(int a, int b) {
            if (lookups[0][a] == lookups[0][b]) {
                return lookups[1][a] == lookups[1][b] ? 0 : (lookups[1][a] < lookups[1][b] ? -1 : 1);
            }
            return lookups[0][a] < lookups[0][b] ? -1 : 1;
        }
    };
}
