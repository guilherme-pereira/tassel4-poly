/*
 *  HDF5ByteGenotype
 */
package net.maizegenetics.pal.alignment.genotype;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import net.maizegenetics.pal.alignment.HapMapHDF5Constants;

/**
 *
 * @author Terry Casstevens
 */
public class HDF5ByteGenotype extends AbstractGenotype {

    private static final int SHIFT_AMOUNT = 16;
    /**
     * Byte representations of DNA sequences are stored in blocks of 65536 sites
     */
    private static final int HDF5_GENOTYPE_BLOCK_SIZE = 1 << SHIFT_AMOUNT;
    private final IHDF5Reader myHDF5Reader;
    private final CacheLoader<Long, byte[]> myGenoLoader = new CacheLoader<Long, byte[]>() {
        public byte[] load(Long key) {
            long offset = getSiteStartFromKey(key) << SHIFT_AMOUNT;
            byte[] data;
            synchronized (myHDF5Reader) {
                data = myHDF5Reader.readAsByteArrayBlockWithOffset(getTaxaGenoPath(getTaxonFromKey(key)), HDF5_GENOTYPE_BLOCK_SIZE, offset);
            }
            return data;
        }
    };
    private final LoadingCache<Long, byte[]> myGenoCache;

    private static long getCacheKey(int taxon, int site) {
        return ((long) taxon << 32) + (site / HDF5_GENOTYPE_BLOCK_SIZE);
    }

    private static int getTaxonFromKey(long key) {
        return (int) (key >>> 32);
    }

    private static int getSiteStartFromKey(long key) {
        return (int) ((key << 32) >>> 32);
    }

    private String getTaxaGenoPath(int taxon) {
        return HapMapHDF5Constants.GENOTYPES + "/taxon" + taxon;
    }

    private HDF5ByteGenotype(IHDF5Reader reader, int numTaxa, int numSites, boolean phased, String[][] alleleEncodings) {
        super(numTaxa, numSites, phased, alleleEncodings);
        myHDF5Reader = reader;
        myGenoCache = CacheBuilder.newBuilder()
                .maximumSize((3 * getTaxaCount()) / 2)
                .build(myGenoLoader);
    }

    static HDF5ByteGenotype getInstance(IHDF5Reader reader) {
        return new HDF5ByteGenotype(reader, 0, 0, false, null);
    }

    @Override
    public byte getBase(int taxon, int site) {
        long key = getCacheKey(taxon, site);
        try {
            byte[] data = myGenoCache.get(key);
            return data[site % HDF5_GENOTYPE_BLOCK_SIZE];
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException("HDF5ByteGenotype: getBase: Error getting base from cache.");
        }
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
