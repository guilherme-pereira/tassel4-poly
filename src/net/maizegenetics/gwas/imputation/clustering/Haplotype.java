package net.maizegenetics.gwas.imputation.clustering;

import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;

/**
 * A Haplotype is a sequence which can be either a true haplotype or a genotype.
 * The natural ordering is by number of non-missing values, descending.
 * @author Peter Bradbury
 */
public class Haplotype implements Comparable<Haplotype>{
	/**
	 * The haplotype (or genotype) sequence
	 */
	public byte[] seq;
	
	/**
	 * The number of non missing values in the sequence
	 */
	public int notMissingCount;
	
	/**
	 * The length of the sequence
	 */
	public int seqlen;
	
	/**
	 * The byte value for missing
	 */
	public static final byte N = NucleotideAlignmentConstants.getNucleotideDiploidByte("N");
	
	/**
	 * The index of the taxon from the originating alignment from which this sequence was taken
	 */
	public int taxonIndex;

	/**
	 * @param hap	the haplotype or genotype sequence
	 * @param taxon	the taxon index
	 */
	public Haplotype(byte[] hap, int taxon) {
		seq = hap;
		seqlen = seq.length;
		countNotMissing();
		taxonIndex = taxon;
	}
	
	/**
	 * @param hap	the haplotype or genotype sequence
	 */
	public Haplotype(byte[] hap) {
		this(hap, -1);
	}
	
	/**
	 * 
	 */
	public void countNotMissing() {
		notMissingCount = 0;
		for (byte b:seq) if (b != N) notMissingCount++;
	}

	@Override
	public int compareTo(Haplotype arg0) {
		if (notMissingCount > arg0.notMissingCount) return -1;
		if (notMissingCount < arg0.notMissingCount) return 1;
		return 0;
	}
	
	/**
	 * Distance is defined as the sum of the distances between each pair of sites. If one or both of the sites are missing then
	 * the distance is zero. If the sites are equal the distance is zero. If both sites are homozygous but different, the distance is 2.
	 * If one site is heterozygous, the distance is 1.
	 * @param h0	another Haplotype. Both haplotypes must have the same sequence length.
	 * @return	the distance between the haplotypes 
	 */
	public int distanceFrom(Haplotype h0) {
		return getDistance(seq, h0.seq);
	}

	/**
	 * @param hap0	the first Haplotype
	 * @param hap1	the second Haplotype
	 * @return distance as defined by distanceFrom(Haplotype h0)
	 */
	public static int getDistance(Haplotype hap0, Haplotype hap1) {
		return hap0.distanceFrom(hap1);
	}
	
	/**
	 * @param hap0	the first Haplotype
	 * @param hap1	the second Haplotype
	 * @return distance as defined by distanceFrom(Haplotype h0)
	 * @see #distanceFrom(Haplotype h0)
	 */
	public static int getDistance(byte[] hap0, byte[] hap1) {
		int d = 0;
		int n = hap0.length;
		for (int s = 0; s < n; s++) {
			byte b0 = hap0[s];
			byte b1 = hap1[s];
			if (b0 != b1 && b0 != N && b1 != N) {
				if (AlignmentUtils.isHeterozygous(b0)) {
					if (!AlignmentUtils.isHeterozygous(b1)) d++;
				} else if (AlignmentUtils.isHeterozygous(b1)) d++;
				else d += 2;
			}
		}
		return d;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (byte b : seq) sb.append(NucleotideAlignmentConstants.getNucleotideIUPAC(b));
		return sb.toString();
	}
	
}
