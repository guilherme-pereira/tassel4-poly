package net.maizegenetics.gwas.imputation.clustering;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;

/**
 * HaplotypeClusterer clusters haplotypes by distance as defined in the Haplotype class.
 * This class can be extended to apply different clustering rules.
 * @author Peter Bradbury
 */

public class HaplotypeClusterer {
	ArrayList<Haplotype> haplotypeList;
	int numberOfHaplotypes;
	int nsites;
	static final byte N = NucleotideAlignmentConstants.getNucleotideDiploidByte("N");
	ArrayList<HaplotypeCluster> clusterList;
	public static enum TYPE {whole, partial}
	
	/**
	 * @param haplotypes	a List of haplotype or genotype sequence to be clustered
	 */
	public HaplotypeClusterer(List<byte[]> haplotypes) {
		numberOfHaplotypes = haplotypes.size();
		nsites = haplotypes.get(0).length;
		haplotypeList = new ArrayList<Haplotype> ();
		for (byte[] hap:haplotypes) haplotypeList.add(new Haplotype(hap));
		Collections.sort(haplotypeList);
	}
	
	/**
	 * @param haplotypes	a List of Haplotypes to be clustered
	 */
	public HaplotypeClusterer(ArrayList<Haplotype> haplotypes) {
		numberOfHaplotypes = haplotypes.size();
		nsites = haplotypes.get(0).seqlen;
		haplotypeList = haplotypes;
		Collections.sort(haplotypeList);
	}

	/**
	 * Groups a list of Haplotypes into clusters.
	 * Clusters are created so that all Haplotypes in a cluster have 0 pairwise distance.
	 * Because of missing data a Haplotype can be assigned to more than one cluster.
	 */
	public void makeClusters() {
		LinkedList<Haplotype> tmpHapList = new LinkedList<Haplotype>(haplotypeList);
		clusterList = new ArrayList<HaplotypeCluster>();
		HaplotypeCluster cluster = new HaplotypeCluster(tmpHapList.removeFirst(), 1.0);
		clusterList.add(cluster);
		
		Iterator<Haplotype> hit = tmpHapList.iterator();
		while (hit.hasNext()) {
			boolean inCluster = false;
			Haplotype hap = hit.next();
			for (HaplotypeCluster clus : clusterList) {
				if (clus.get(0).distanceFrom(hap) == 0) {
					inCluster = true;
					break;
				}
			}
			if (!inCluster) {
				HaplotypeCluster newCluster = new HaplotypeCluster(hap, 1.0);
				clusterList.add(newCluster);
				hit.remove();
			}
		}
		
		int nclusters = clusterList.size();
		for (Haplotype hap:tmpHapList) {
			boolean[] incluster = new boolean[nclusters];
			Arrays.fill(incluster,true);
			int count = 0;
			for (int c = 0; c < nclusters; c++) {
				HaplotypeCluster clus = clusterList.get(c);
				Iterator<Haplotype> clusIt = clus.getIterator();
				while(clusIt.hasNext() && incluster[c]) {
					Haplotype member = clusIt.next();
					if (member.distanceFrom(hap) > 0) {
						incluster[c] = false;
					}
				}
				if (incluster[c]) {
					count++;
					clus.add(hap);
				}
			}
			double hapscore = 1.0 / ((double) count);
			for (int c = 0; c < nclusters; c++) {
				if (incluster[c]) clusterList.get(c).incrementScore(hapscore);
			}
		}
	}
	
	/**
	 * Removes all haplotypes within maxdistance of the haplotype of the first cluster.
	 * After the haplotypes have been removed clusters are remade and sorted.
	 * @param maxdistance	all haplotypes within maxdistance or less of the haplotype of the first cluster are removed from the haplotypeList
	 */
	public HaplotypeCluster removeFirstHaplotypes(int maxdistance) {
		
		//remove all haplotypes maxdistance or less from the haplotype of cluster 0
		Iterator<Haplotype> hapit = haplotypeList.listIterator();
		byte[] firstHaplotype = clusterList.get(0).getHaplotype();
		ArrayList<Haplotype> haplist = new ArrayList<Haplotype>();
		while (hapit.hasNext()) {
			Haplotype testhap = hapit.next();
			if (Haplotype.getDistance(firstHaplotype, testhap.seq) <= maxdistance) {
				haplist.add(testhap);
				hapit.remove();
			}
		}

		makeClusters();
		sortClusters();
		return new HaplotypeCluster(haplist, haplist.size());
	}
	
	/**
	 * Sorts clusters according to HaplotypeCluster sort order.
	 */
	public void sortClusters() {
		Collections.sort(clusterList);
	}
	
	/**
	 * @return	the number of clusters
	 */
	public int getNumberOfClusters() {
		return clusterList.size();
	}
	
	/**
	 * @return the number of Haplotypes contained in each cluster
	 */
	public int[] getClusterSizes() {
		int nclusters = clusterList.size();
		int[] sizes = new int[nclusters];
		for (int i = 0; i < nclusters; i++) sizes[i] = clusterList.get(i).getSize();
		return sizes;
	}
	
	/**
	 * After the initial cluster formation a Haplotype score equals the 1 / (number of clusters to which it belongs).
	 * Merging does not update the cluster score.
	 * @return	the score for each cluster.
	 */
	public double[] getClusterScores() {
		int n = clusterList.size();
		double[] scores = new double[n];
		for (int i = 0; i < n; i++) {
			scores[i] = clusterList.get(i).getScore();
		}
		return scores;
	}
	
	/**
	 * @return	List of clusters
	 */
	public ArrayList<HaplotypeCluster> getClusterList() { return clusterList; }
	
	/**
	 * @param cluster0	
	 * @param cluster1	
	 * @return the minimum of the number of haplotypes in cluster0 not in cluster1 and the number of haplotypes in cluster 1 not in cluster0
	 */
	public static int clusterDistanceDistinctHaplotypes(HaplotypeCluster cluster0, HaplotypeCluster cluster1){
		int count0 = cluster0.getCountOfHaplotypesNotInThisCluster(cluster1);
		int count1 = cluster1.getCountOfHaplotypesNotInThisCluster(cluster0);
		return Math.min(count0, count1);
	}
	
	/**
	 * @param cluster0	
	 * @param cluster1	
	 * @return the minimum of the number of haplotypes in cluster0 not in cluster1 and the number of haplotypes in cluster 1 not in cluster0 divided my the minimum number of haplotypes in the two clusters
	 */
	public static double clusterDistanceDistinctHaplotypeProportion(HaplotypeCluster cluster0, HaplotypeCluster cluster1){
		int count0 = cluster0.getCountOfHaplotypesNotInThisCluster(cluster1);
		int count1 = cluster1.getCountOfHaplotypesNotInThisCluster(cluster0);
		double minNumber = Math.min(cluster0.getSize(), cluster1.getSize());
		return Math.min(count0, count1) / minNumber;
	}
	
	/**
	 * @param cluster0
	 * @param cluster1
	 * @return the number of alleles difference between the cluster haplotypes
	 */
	public static int clusterDistanceClusterHaplotypeDiff(HaplotypeCluster cluster0, HaplotypeCluster cluster1){
		byte[] hap0 = cluster0.getHaplotype();
		byte[] hap1 = cluster1.getHaplotype();
		return Haplotype.getDistance(hap0, hap1);
	}
	
	/**
	 * @param cluster0
	 * @param cluster1
	 * @return the number of alleles difference between the cluster haplotypes divided by the number of nonmissing sites
	 */
	public static double clusterDistanceClusterDiffProportion(HaplotypeCluster cluster0, HaplotypeCluster cluster1){
		byte[] hap0 = cluster0.getHaplotype();
		byte[] hap1 = cluster1.getHaplotype();
		int n = hap0.length;
		double notmissing = 0;
		double notmatching = 0;
		byte N = NucleotideAlignmentConstants.getNucleotideDiploidByte('N');
		for (int i = 0; i < n; i++) {
			if (hap0[i] != N && hap1[i] != N) {
				notmissing++;
				if (hap0[i] != hap1[i]) notmatching++;
			}
		}
		return notmatching/notmissing;
	}
	
	/**
	 * @param cluster0
	 * @param cluster1
	 * @return the maximum of the pairwise difference between individual haplotypes in the two clusters
	 */
	public static int clusterDistanceMaxPairDiff(HaplotypeCluster cluster0, HaplotypeCluster cluster1){
		HaplotypeCluster subCluster0 = cluster0.copy();
		subCluster0.removeAll(cluster1);;
		HaplotypeCluster subCluster1 = cluster1.copy();
		subCluster1.removeAll(cluster0);
		
		int maxdiff = 0;
		
		Iterator<Haplotype> hit0 = subCluster0.getIterator();
		while (hit0.hasNext()) {
			Haplotype h0 = hit0.next();
			Iterator<Haplotype> hit1 = subCluster1.getIterator();
			while (hit1.hasNext()) {
				Haplotype h1 = hit1.next();
				maxdiff = Math.max(maxdiff, Haplotype.getDistance(h0, h1));
			}
		}
		return maxdiff;
	}
	
	/**
	 * @param cluster0
	 * @param cluster1
	 * @return the average of the pairwise differences between the haplotypes that are not shared between clusters
	 */
	public static double clusterDistanceAveragePairDiff(HaplotypeCluster cluster0, HaplotypeCluster cluster1){
		HaplotypeCluster subCluster0 = cluster0.copy();
		subCluster0.removeAll(cluster1);;
		HaplotypeCluster subCluster1 = cluster1.copy();
		subCluster1.removeAll(cluster0);
		
		int totaldiff = 0;
		int count = 0;
		
		Iterator<Haplotype> hit0 = subCluster0.getIterator();
		while (hit0.hasNext()) {
			Haplotype h0 = hit0.next();
			Iterator<Haplotype> hit1 = subCluster1.getIterator();
			while (hit1.hasNext()) {
				Haplotype h1 = hit1.next();
				totaldiff += Haplotype.getDistance(h0, h1);
				count++;
			}
		}
		return ((double) totaldiff) / ((double) count);
	}
	
	/**
	 * @param cluster0
	 * @param cluster1
	 * @return	the total of the pairwise differences between the haplotypes that are not shared between clusters
	 */
	public static int clusterDistanceTotalPairDiff(HaplotypeCluster cluster0, HaplotypeCluster cluster1){
		HaplotypeCluster subCluster0 = cluster0.copy();
		subCluster0.removeAll(cluster1);;
		HaplotypeCluster subCluster1 = cluster1.copy();
		subCluster1.removeAll(cluster0);
		
		int totaldiff = 0;
		
		Iterator<Haplotype> hit0 = subCluster0.getIterator();
		while (hit0.hasNext()) {
			Haplotype h0 = hit0.next();
			Iterator<Haplotype> hit1 = subCluster1.getIterator();
			while (hit1.hasNext()) {
				Haplotype h1 = hit1.next();
				totaldiff += Haplotype.getDistance(h0, h1);
			}
		}
		return totaldiff;
	}
	
	/**
	 * Merges clusters whose maximum pairwise difference is less than maxdiff. Clusters are tested sequentially.
	 * That is, if two clusters are merged, they become the new head cluster against which remaining clusters are tested for merging.
	 * @param candidateClusters	an ArrayList of HaplotypeClusters
	 * @param maxdiff
	 * @return merges clusters whose maximum pairwise difference is less or equal to maxdiff
	 */
	public static ArrayList<HaplotypeCluster> getMergedClusters(ArrayList<HaplotypeCluster> candidateClusters, int maxdiff) {
		ArrayList<HaplotypeCluster> candidates = new ArrayList<HaplotypeCluster>(candidateClusters);
		Collections.sort(candidates);
		ArrayList<HaplotypeCluster> mergedClusterList = new ArrayList<HaplotypeCluster>();

		while (candidates.size() > 0) {
			HaplotypeCluster headCluster = candidates.remove(0);
			mergedClusterList.add(headCluster);
			Iterator<HaplotypeCluster> cit = candidates.iterator();
			while (cit.hasNext()) {
				HaplotypeCluster candidate = cit.next();
				boolean mergePair = doMerge(headCluster, candidate, maxdiff);
				if (mergePair) {
					cit.remove();
					mergeTwoClusters(headCluster, candidate);
				}
			} 
		}
		return mergedClusterList;
	}

	/**
	 * Tests whether two clusters are less than or equal to maxdiff distant. 
	 * Uses clusterDistanceMaxPairDiff() to calculate distance.
	 * @param c0	a cluster
	 * @param c1	another cluster
	 * @param maxdiff	
	 * @return	true if clusters are maxdiff or less distant, false otherwise
	 */
	public static boolean doMerge(HaplotypeCluster c0, HaplotypeCluster c1, int maxdiff) {
		int diff = clusterDistanceMaxPairDiff(c0, c1);
		return diff <= maxdiff;
	}
	
	/**
	 * Merges two clusters, c0 and c1.
	 * @param c0
	 * @param c1
	 */
	public static void mergeTwoClusters(HaplotypeCluster c0, HaplotypeCluster c1) {
		c1.removeAll(c0);
		c0.addAll(c1);
	}
}
