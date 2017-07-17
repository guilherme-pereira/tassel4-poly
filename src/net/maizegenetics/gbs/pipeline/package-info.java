/**
 * The GBS package provide plugins and analysis approaches for genotyping by sequencing.
 * <p>
 * Genotyping-by-sequencing (GBS) is a reduced representation approach for sampling the genome.
 * The approach dates back to 2000 and 2001, when the human genome was restricted
 * and sequenced by Sanger.
 * <p>
 * The TASSEL GBS pipeline grew out of a protocol developed out of efficient and 
 * inexpensive approach published in 2011.  
 * The TASSEL GBS PE(pair end) pipeline is used to improve physical alignment of GBS tags.
 * <p>
 * Elshire RJ, Glaubitz JC, Sun Q, Poland JA, Kawamoto K, Buckler ES, Mitchell SE. 
 * (2011) A robust, simple genotyping-by-sequencing (GBS) approach for high diversity species. 
 * PLoS One 6(5): e19379.
 * 
 * The bioinformatics has developed substantially since then and a new publication describing
 * the approaches is forthcoming:
 * <p>
 * Glaubitz, Casstevens, Harriman, Elshire, and Buckler (2013) in prep.
 * <p>
 * Key design principles:
 * Observed sequences are calls READS.  These reads are trimmed and combined into
 * clusters of identical sequences called TAGS.
 * The distribution of tags are scored across taxa (samples) and recorded in a Tags by Taxa object (TBT).
 * The genetic and physical mapping of tags are recorded in a Tags On Physical Map (TOPM)
 * object.  Once SNPs are identified their positions are recorded in a TOPM, which 
 * can be used in a single step production calling pipeline.
 *
 * <p>
 * This effort was funded by the NSF Plant Genome, NSF BREAD program, and the USDA.
 * 
 */
package net.maizegenetics.gbs.pipeline;
