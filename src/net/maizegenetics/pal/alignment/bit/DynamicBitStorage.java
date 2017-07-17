package net.maizegenetics.pal.alignment.bit;

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import net.maizegenetics.pal.alignment.AlignmentNew.ALLELE_SCOPE_TYPE;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.genotype.Genotype;
import net.maizegenetics.util.BitSet;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutionException;

/**
 * Provides rapid conversion routines and caching from byte encoding of
 * nucleotides to bit encoding. Only two alleles are supported for each scope
 * (e.g. Major & minor, or Reference & Alternate).
 * <p></p>
 * The cache is designed to support multiple scopes, but currently scope must be
 * passed in at construction.
 * <p></p>
 * It is not clear that site or taxa optimization is needed. The code should be
 * highly parallelizable as long as the gets are not for adjacent sites.
 *
 * @author Ed Buckler
 */
public class DynamicBitStorage implements BitStorage {

    private Genotype myGenotype;
    private ALLELE_SCOPE_TYPE myPreferredScope = ALLELE_SCOPE_TYPE.Frequency;
    private byte[] myPrefAllele0;  //usually 0 is major or reference
    private byte[] myPrefAllele1;  //usually 1 is minor or alternate
    private final int myTaxaCount;
    private final int mySiteCount;
    public static final int SBoff = 58;
    public static final int AllScopeoff = 52;

    private enum SB {

        TAXA(0), SITE(1);
        public final int index;

        SB(int index) {
            this.index = index;
        }
    };
    private LoadingCache<Long, BitSet[]> bitCache; //taxon id number, BitSet[2] = {MajorBitSet, MinorBitSet}
    private CacheLoader<Long, BitSet[]> bitLoader = new CacheLoader<Long, BitSet[]>() {
        public BitSet[] load(Long key) {
            BitSet[] bs;
            if (getDirectionFromKey(key) == SB.TAXA) {
                byte[] a1 = myPrefAllele0;
                byte[] a2 = myPrefAllele1;
                int taxon = getSiteOrTaxonFromKey(key);
                bs = AlignmentUtils.calcBitPresenceFromGenotype(myGenotype.getBaseRow(taxon), a1, a2); //allele comp
                return bs;
            } else {
                ArrayList toFill = new ArrayList<Integer>();
                toFill.add(key);
                try {
                    bitCache.putAll(loadAll(toFill));
                    return bitCache.get(key);
                } catch (Exception e) {
                    e.printStackTrace();
                    return null;
                }
            }
        }

        @Override
        public Map<Long, BitSet[]> loadAll(Iterable<? extends Long> keys) throws Exception {
            long key = keys.iterator().next();
            //This pivoting code is needed if myGenotype is store in taxa direction
            //It runs about 7 times faster than getting base sequentially across taxa.
            HashMap<Long, BitSet[]> result = new HashMap<Long, BitSet[]>(64);
            int site = getSiteOrTaxonFromKey(key);
            int length = (mySiteCount - site < 64) ? mySiteCount - site : 64;
            byte[][] genotypeTBlock = new byte[length][myTaxaCount];
            for (int t = 0; t < myTaxaCount; t++) {
                for (int s = 0; s < genotypeTBlock.length; s++) {
                    genotypeTBlock[s][t] = myGenotype.getBase(t, site + s);
                }
            }
            for (int i = 0; i < length; i++) {
                byte a1 = myPrefAllele0[site + i];
                byte a2 = myPrefAllele1[site + i];
                BitSet[] bs = AlignmentUtils.calcBitPresenceFromGenotype(genotypeTBlock[i], a1, a2);
                result.put(getKey(SB.SITE, myPreferredScope, site + i), bs);
            }
            return result;
        }
    };

    private long getKey(SB direction, ALLELE_SCOPE_TYPE aT, int siteOrTaxon) {
        return ((long) direction.index << SBoff) | ((long) aT.ordinal() << AllScopeoff) | (long) siteOrTaxon;
    }

    private int getSiteOrTaxonFromKey(long key) {
        return (int) ((key << 32) >>> 32);
    }

    private SB getDirectionFromKey(long key) {
        if (key >>> SBoff == SB.TAXA.index) {
            return SB.TAXA;
        }
        return SB.SITE;
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        try {
            return bitCache.get(getKey(SB.TAXA, myPreferredScope, taxon))[alleleNumber];
        } catch (ExecutionException e) {
            e.printStackTrace();
            return null;
        }
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        try {
            return bitCache.get(getKey(SB.SITE, myPreferredScope, site))[alleleNumber];
        } catch (ExecutionException e) {
            e.printStackTrace();
            return null;
        }
    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        BitSet result = getAllelePresenceForAllSites(taxon, alleleNumber);
        if (result == null) {
            return new long[0];
        }
        return result.getBits(startBlock, endBlock - 1);   //BitSet is inclusive, while this method is exclusive.
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber) {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber) {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    @Override
    public long[] getPhasedAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock) {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    public DynamicBitStorage(Genotype genotype, ALLELE_SCOPE_TYPE currentScope, byte[] prefAllele0, byte[] prefAllele1) {
        myGenotype = genotype;
        myPreferredScope = currentScope;
        mySiteCount = myGenotype.getSiteCount();
        myTaxaCount = myGenotype.getTaxaCount();
        myPrefAllele0 = Arrays.copyOf(prefAllele0, prefAllele0.length);
        myPrefAllele1 = Arrays.copyOf(prefAllele1, prefAllele1.length);
        bitCache = CacheBuilder.newBuilder()
                .maximumSize(3_000_000)
                .build(bitLoader);
    }

    public static DynamicBitStorage getInstance(Genotype genotype, ALLELE_SCOPE_TYPE currentScope, byte[] prefAllele) {
        int numSites = prefAllele.length;
        byte[] prefAllele0 = new byte[numSites];
        byte[] prefAllele1 = new byte[numSites];
        for (int i = 0; i < numSites; i++) {
            prefAllele0[i] = (byte) (prefAllele[i] >>> 4);
            prefAllele1[i] = (byte) (prefAllele[i] & 0xf);
        }
        return new DynamicBitStorage(genotype, currentScope, prefAllele0, prefAllele1);
    }
}
