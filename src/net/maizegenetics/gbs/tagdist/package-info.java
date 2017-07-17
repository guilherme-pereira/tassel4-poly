/**
 * Basic data structures for holding tags.
 *<p>
 * Definitions:
 * <li>Read:  a single sequence from a sequencer
 * <li>Tag:  a set of reads with a unique sequence
 * <li>Taxa:  a sample or operational taxonomic unit
 * <li>Counts:  the number of reads seen for a particular tag
 * <li>PE Tag:  a pair of pair end tags, including the forward tag and the backward tag
 * <li>Contig:  a contig is formed by a PE Tag whose forward tag and backward tag have overlap 
 *<p>
 * All sequences are encoded using BaseEncoder.java, which compresses each base into
 * two bit representations.  Tags are also encoded in generally the length of two
 * longs (64bits x 2 = 128 bits / 2bitsPerBase = 64 bases).
 * PE Tag has 8 longs (32 * 8 = 256 bases)
 *<p>
 * If tags are shorter than current length, then they are padded with polyA.
 * The length of the tag is recorded in TagLength.
 *<p>
 * All the basic data structures inherit from {@link net.maizegenetics.gbs.tagdist.Tags}.
 * {@link net.maizegenetics.gbs.tagdist.TagCounts} adds
 * information on the number of reads.  
 * {@link net.maizegenetics.gbs.tagdist.TagsByTaxa} add information on the
 * distributions of Tags across taxa.
 *
 */
package net.maizegenetics.gbs.tagdist;
