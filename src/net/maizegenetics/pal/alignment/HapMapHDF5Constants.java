/*
 * HapMapHDF5Constants
 */
package net.maizegenetics.pal.alignment;

/**
 *
 * @author terry
 */
public final class HapMapHDF5Constants {
    // Paths
    public static final String ROOT = "/";
    public static final String TAXA = "Taxa";
    public static final String ALLELE_STATES = "AlleleStates";
    public static final String POSITIONS = "Positions";
    public static final String ALLELES = "Alleles";  //TODO - like to move this into SITE_DESC

    public static final String TBIT = "TBit";  //Taxa optimized bit alignment
    public static final String SBIT = "SBit";  //Site optimized bit alignment
    public static final String DEPTH = "Depth";  //taxa optimized depth of the base calls
    public static final String GENOTYPES = "Genotypes";  //Taxa based byte alignment
    public static final String LOCI = "SeqRegion";
    public static final String LOCUS_OFFSETS = "SeqRegionOffsets";
    public static final String LOCUS_INDICES = "SeqRegionIndices";
    public static final String SNP_IDS = "SnpIds";
    public static final String SITE_DESC = "SiteDesc/";
 
    public static final String ALLELE_CNT = SITE_DESC+"AlleleCnt";
    public static final String MAF = SITE_DESC+"MAF";
    public static final String SITECOV = SITE_DESC+"SiteCoverage";
    public static final String SITEHET = SITE_DESC+"SiteHET";
 //   public static final String MINOR_ALLELE_FREQUENCY = "MAF";
    public static final String MAJOR_ALLELE = SITE_DESC+"MajorAllele";
    public static final String MINOR_ALLELE = SITE_DESC+"MinorAllele";
    public static final String ALLELE_FREQ_ORD = SITE_DESC+"AlleleFreqOrder";
    public static final String REF_ALLELE = SITE_DESC+"REFAllele";
    
    public static final String TAXA_DESC = "TaxaDesc/";
    public static final String TAXACOV = TAXA_DESC+"TaxaCoverage";
    public static final String TAXAHET = TAXA_DESC+"TaxaHet";
    
    public static final String LD_DESC = SITE_DESC+"/LD";
    public static final String LDR2_DESC = LD_DESC+"/R2";
    public static final String LDP_DESC = LD_DESC+"/P";
    public static final String LDPropLD_DESC = LD_DESC+"/PropLD";
    public static final String LDMinDist_DESC = LD_DESC+"/MinDist";
    public static final String BPLDMean34 = LD_DESC+"/BPLDMean34R2";
    public static final String BPLDMeanMax = LD_DESC+"/BPLDMeanMaxR2";
    public static final String BPmaxMaxLD = LD_DESC+"/BPmaxMaxLDR2";
    public static final String BPminMaxLD = LD_DESC+"/BPminMaxLDR2";
    public static final String BPPopCnt = LD_DESC+"/BPPopCnt";
    
    
    public static final String ERROR_DESC = SITE_DESC+"/ERROR";
    public static final String BPECERROR_DESC = ERROR_DESC+"/BPECERROR";
    public static final String BPECAVGR2_DESC = ERROR_DESC+"/BPECR2";
    
    public static final String IBSERRORRATE_DESC = ERROR_DESC+"/IBSErrorRate";
    public static final String IBSMINORERRORRATE_DESC = ERROR_DESC+"/IBSMinorErrorRate";
    public static final String IBSERROR_DESC = ERROR_DESC+"/IBSError";
    public static final String IBSMINORCORR_DESC = ERROR_DESC+"/IBSMinorCorrect";
    public static final String IBSMAJORCORR_DESC = ERROR_DESC+"/IBSMajorCorrect";
    
    // Attributes
    public static final String DEFAULT_ATTRIBUTES_PATH = ROOT;
    public static final String NUM_SBIT_WORDS = "numSBitWords";
    public static final String NUM_TBIT_WORDS = "numTBitWords";
    public static final String NUM_TAXA = "numTaxa";
    public static final String NUM_SITES = "numSites";
    public static final String MAX_NUM_ALLELES = "maxNumAlleles";
    public static final String RETAIN_RARE_ALLELES = "retainRareAlleles";
    public static final String NUM_LD_BINS = LD_DESC+"/numLDBins";
    public static final String LD_BINS = LD_DESC+"/binsLD";



    // Annotations
    public static final String GWAS = "GWAS";
    public static final String GENOMIC = "GenomeAnno";
    public static final String POP_GEN = "PopGenAnno";


    
    private HapMapHDF5Constants() {
        // do not instantiate
    }
    
}