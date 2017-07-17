/*
 *  AlleleDepth
 */
package net.maizegenetics.pal.alignment.depth;

/**
 *
 * @author terry, Gabriel Rodrigues Alves Margarido
 */
public interface AlleleDepth {

    /**
     * Returns depth count for each diploid allele at the given taxon and site.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return two counts
     */
    public short[] getDepthForAlleles(int taxon, int site);
}
