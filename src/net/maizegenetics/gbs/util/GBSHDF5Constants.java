
package net.maizegenetics.gbs.util;

/**
 * Path definitions for the HDF5 GBS file formats.
 * 
 * @author edbuckler
 */
public class GBSHDF5Constants {
    // Paths
    public static final String ROOT = "/";
    public static final String TAGS="tags";
    public static final String TAGLENGTH="tagLength";
    public static final String MULTIMAPS="multimaps";
    public static final String BEST_STRAND="bestStrand";
    public static final String BEST_CHR="bestChr";
    public static final String BEST_STARTPOS="bestStartPos";
    public static final String BEST_ENDPOS="bestEndPos";
    public static final String BEST_DIVERGENCE="bestDivergence";
    public static final String BEST_MAPP="bestMapP";
    public static final String BEST_DCOP="bestDcoP";
    public static final String BEST_EVIDENCE="bestEvidence";
    public static final String BEST_MAPINDEX="bestMapIndex";
    public static final String GENETICMAMMPING="gm";
    public static final String GENETICMAMMPINGGW="gmgw";
    public static final String MAPBASE="map";
    public static final String VARIANTDEF="variantDef";
    public static final String VARIANTPOSOFF="variantPosOff";
    
    //Attributes
    public static final String TAGCOUNT = "tagCount";
    public static final String MAXVARIANTS = "maxVariants";
    public static final String MAXMAPPING = "maxMapping";
    public static final String TAGLENGTHINLONG = "tagLengthInLong";
    public static final String IFHASGM = "ifHasGM";
    public static final String IFHASGMGW = "ifHasGMGW";
    
    private GBSHDF5Constants() {
    }
    
    
    
}
