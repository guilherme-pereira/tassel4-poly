package net.maizegenetics.gwas.imputation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Random;
import java.util.Set;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

import net.maizegenetics.baseplugins.ConvertSBitTBitPlugin;
import net.maizegenetics.gwas.NAM.AGPMap;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.util.BitSet;

public class ImputationUtils {
	private static final Logger myLogger = Logger.getLogger(ImputationUtils.class);

	public static Pattern tab = Pattern.compile("\t");
	
	public static int[] order(int[] array) {
		class SortElement implements Comparable<SortElement> {
			int val;
			int ndx;
			
			SortElement(int x, int index) {
				val = x;
				ndx = index;
			}
			
			@Override
			public int compareTo(SortElement se) {
				return val - se.val;
			}
		}

		int n = array.length;
		SortElement[] sortArray = new SortElement[n];
		for (int i = 0; i < n; i++) {
			sortArray[i] = new SortElement(array[i], i);
		}
		
		Arrays.sort(sortArray);
		int[] order = new int[n];
		for (int i = 0; i < n; i++) {
			order[i] = sortArray[i].ndx;
		}
		return order;
	}
	
	public static int[] reverseOrder(int[] array) {
		class SortElement implements Comparable<SortElement> {
			int val;
			int ndx;
			
			SortElement(int x, int index) {
				val = x;
				ndx = index;
			}
			
			@Override
			public int compareTo(SortElement se) {
				return se.val - val;
			}
		}

		int n = array.length;
		SortElement[] sortArray = new SortElement[n];
		for (int i = 0; i < n; i++) {
			sortArray[i] = new SortElement(array[i], i);
		}
		
		Arrays.sort(sortArray);
		int[] order = new int[n];
		for (int i = 0; i < n; i++) {
			order[i] = sortArray[i].ndx;
		}
		return order;
	}
	
	public static Alignment[] getTwoClusters(Alignment a, int[] parentIndex) {
		int maxiter = 5;
		Alignment tb = ConvertSBitTBitPlugin.convertAlignment(a, ConvertSBitTBitPlugin.CONVERT_TYPE.tbit, null);
		//if (a instanceof TBitAlignment) tb = (TBitAlignment) a;
		//else tb = TBitAlignment.getInstance(a);
		
		//if the parents are in the data set use these as seeds
		//if one parent is in the dataset pick the taxon farthest from it as the other seed
		//if neither parent is in the dataset choose random seeds
		int ntaxa = tb.getSequenceCount();
		int nsnps = tb.getSiteCount();
		int seed1 = parentIndex[0];
		int seed2 = parentIndex[1];
		float[] loc1, loc2;
		Random rand = new Random();
		if (seed1 == -1) {
			if (seed2 == -1) {
				//both parents are not in the data set
				seed1 = rand.nextInt(ntaxa);
				loc1 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(seed1, 0), tb.getAllelePresenceForAllSites(seed1, 1)}, nsnps);
				while (seed2 == -1 || seed1 == seed2) {
					seed2 = rand.nextInt(ntaxa);
				}
				loc2 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(seed2, 0), tb.getAllelePresenceForAllSites(seed2, 1)}, nsnps);
			} else {
				//parent2 is in the data set
				loc1 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(0, 0), tb.getAllelePresenceForAllSites(0, 1)}, nsnps);
				loc2 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(seed2, 0), tb.getAllelePresenceForAllSites(seed2, 1)}, nsnps);
				seed1 = 0;
				float prevdist = getManhattanDistance(loc2, loc1, nsnps);
				for (int t = 1; t < ntaxa; t++) {
					float[] tloc = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(t, 0), tb.getAllelePresenceForAllSites(t, 1)}, nsnps);
					float dist = getManhattanDistance(loc2, tloc, nsnps);
					if (dist > prevdist) {
						prevdist = dist;
						loc1 = tloc;
						seed1 = t;
					}
				}
			}
		} else if (seed2 == -1) {
			//parent1 is in the data set
			loc1 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(seed1, 0), tb.getAllelePresenceForAllSites(seed1, 1)}, nsnps);
			loc2 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(0, 0), tb.getAllelePresenceForAllSites(0, 1)}, nsnps);
			seed2 = 0;
			float prevdist = getManhattanDistance(loc1, loc2, nsnps);
			for (int t = 1; t < ntaxa; t++) {
				float[] tloc = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(t, 0), tb.getAllelePresenceForAllSites(t, 1)}, nsnps);
				float dist = getManhattanDistance(loc1, tloc, nsnps);
				if (dist > prevdist) {
					prevdist = dist;
					loc2 = tloc;
					seed2 = t;
				}
			}
		} else {
			//both parents are in the data set
			loc1 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(seed1, 0), tb.getAllelePresenceForAllSites(seed1, 1)}, nsnps);
			loc2 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(seed2, 0), tb.getAllelePresenceForAllSites(seed2, 1)}, nsnps);
		}
		
		int[] size1 = new int[nsnps];
		int[] size2 = new int[nsnps];
		for (int i = 0; i < nsnps; i++) {
			if (loc1[i] >= 0) size1[i] = 1;
			if (loc2[i] >= 0) size2[i] = 1;
		}
		boolean[] isInCluster1 = new boolean[ntaxa];
		isInCluster1[seed1] = true;
		isInCluster1[seed2] = false;
		
		
		//do initial cluster assignment
		for (int t = 0; t < ntaxa; t++) {
			if (t != seed1 && t != seed2) {
				float[] tloc = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(t, 0), tb.getAllelePresenceForAllSites(t, 1)}, nsnps);
				float dist1 = getManhattanDistance(loc1, tloc, nsnps);
				float dist2 = getManhattanDistance(loc2, tloc, nsnps);
				if (dist1 <= dist2) {
					isInCluster1[t] = true;
					loc1 = getMeanLocation(loc1, size1, tloc, true, nsnps);
				} else {
					isInCluster1[t] = false;
					loc2 = getMeanLocation(loc2, size2, tloc, true, nsnps);
				}
			}
		}
		
		//update cluster membership until there are no changes or for the maximum number of iterations
		for (int iter = 0; iter < maxiter; iter++) {
			boolean noChanges = true;
			for (int t = 0; t < ntaxa; t++) {
				float[] tloc = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(t, 0), tb.getAllelePresenceForAllSites(t, 1)}, nsnps);
				float dist1 = getManhattanDistance(loc1, tloc, nsnps);
				float dist2 = getManhattanDistance(loc2, tloc, nsnps);
				if (dist1 <= dist2 && isInCluster1[t] == false) {
					isInCluster1[t] = true;
					loc1 = getMeanLocation(loc1, size1, tloc, true, nsnps);
					loc2 = getMeanLocation(loc2, size2, tloc, false, nsnps);
					noChanges = false;
				} else if (dist1 > dist2 && isInCluster1[t] == true){
					isInCluster1[t] = false;
					loc1 = getMeanLocation(loc1, size1, tloc, false, nsnps);
					loc2 = getMeanLocation(loc2, size2, tloc, true, nsnps);
					noChanges = false;
				}
			}

			if (noChanges) break;
		}
		
		System.out.println("distance between clusters = " + getManhattanDistance(loc1, loc2, nsnps));
		
		//make alignments based on the clusters
		boolean[] isInCluster2 = new boolean[ntaxa];
		for (int t = 0; t < ntaxa; t++) isInCluster2[t] = !isInCluster1[t];
		IdGroup id1 = IdGroupUtils.idGroupSubset(tb.getIdGroup(), isInCluster1);
		IdGroup id2 = IdGroupUtils.idGroupSubset(tb.getIdGroup(), isInCluster2);
		
		Alignment a1 = FilterAlignment.getInstance(tb, id1);
		Alignment a2 = FilterAlignment.getInstance(tb, id2);
		
		return new Alignment[]{a1, a2};
	}
	
	public static Alignment[] getTwoClusters(Alignment inputAlignment, int minGametesPerTaxon) {
		
		//filter out low coverage taxa
		int ntaxa = inputAlignment.getSequenceCount();
		boolean[] include = new boolean[ntaxa];
		
		int nIncluded = 0;
		for (int t = 0; t < ntaxa; t++) {
			if (inputAlignment.getTotalGametesNotMissingForTaxon(t) >= minGametesPerTaxon) {
				nIncluded++;
				include[t] = true;
			}
			else include[t] = false;
		}
		
		Alignment myAlignment;
		if (nIncluded < 10) {
			myLogger.info("Included lines less than 10 in getTwoClusters, poor coverage in interval starting at " + inputAlignment.getSNPID(0));
			return null;
		} else {
			Alignment fa = FilterAlignment.getInstance(inputAlignment, IdGroupUtils.idGroupSubset(inputAlignment.getIdGroup(), include));
			
			myAlignment = BitAlignment.getInstance(fa, false);
		}
		int ntrials = 5;
		int maxiter = 5;
		
		//if the parents are in the data set use these as seeds
		//if one parent is in the dataset pick the taxon farthest from it as the other seed
		//if neither parent is in the dataset choose random seeds
		ntaxa = myAlignment.getSequenceCount();
		int nsnps = myAlignment.getSiteCount();
		boolean[][] isInCluster1 = new boolean[ntrials][ntaxa];
		int bestTrial = -1;
		float maxDistance = 0;
		
		Random rand = new Random();
		
		float[][] taxaLocs = new float[ntaxa][nsnps];
		
		myAlignment.optimizeForTaxa(null);
		for (int t = 0; t < ntaxa; t++) {
			taxaLocs[t] = snpsAsFloatVector(new BitSet[]{myAlignment.getAllelePresenceForAllSites(t, 0), myAlignment.getAllelePresenceForAllSites(t, 1)}, nsnps);
		}
		
		for (int trial = 0; trial < ntrials; trial++) {
			int seed1 = rand.nextInt(ntaxa);
			int seed2 = -1;
			while (seed2 == -1 || seed1 == seed2) {
				seed2 = rand.nextInt(ntaxa);
			}

			isInCluster1[trial][seed1] = true;
			isInCluster1[trial][seed2] = false;
			
			
			//do initial cluster assignment
			for (int t = 0; t < ntaxa; t++) {
				if (t != seed1 && t != seed2) {
					float dist1 = getManhattanDistance(taxaLocs[seed1], taxaLocs[t], nsnps);
					float dist2 = getManhattanDistance(taxaLocs[seed2], taxaLocs[t], nsnps);
					if (dist1 < dist2) {
						isInCluster1[trial][t] = true;
					} else if (dist1 > dist2){
						isInCluster1[trial][t] = false;
					} else if (rand.nextDouble() > 0.5) {
						isInCluster1[trial][t] = true;
					} else {
						isInCluster1[trial][t] = false;
					}
				}
			}
			
			//update cluster membership until there are no changes or for the maximum number of iterations
			float[][] meanLocs = new float[2][];
			boolean badclusters = false;
			for (int iter = 0; iter < maxiter; iter++) {
				boolean noChanges = true;
				
				int nCluster1 = 0;
				int nCluster2 = 0;
				for (int t = 0; t < ntaxa; t++) {
					if (isInCluster1[trial][t]) nCluster1++;
					else nCluster2++;
				}
				
				if (nCluster1 == 0 || nCluster2 == 0) {
					badclusters = true;
					break;
				}
				
				float[][] cluster1Locs = new float[nCluster1][];
				float[][] cluster2Locs = new float[nCluster2][];
				int countCluster1 = 0;
				int countCluster2 = 0;
				for (int t = 0; t < ntaxa; t++) {
					if (isInCluster1[trial][t]) cluster1Locs[countCluster1++] = taxaLocs[t]; 
					else cluster2Locs[countCluster2++] = taxaLocs[t];
				}
				
				meanLocs = new float[2][];
				meanLocs[0] = getMeanLocation(cluster1Locs);
				meanLocs[1] = getMeanLocation(cluster2Locs);
				for (int t = 0; t < ntaxa; t++) {
					float[] tloc = snpsAsFloatVector(new BitSet[]{myAlignment.getAllelePresenceForAllSites(t, 0), myAlignment.getAllelePresenceForAllSites(t, 1)}, nsnps);
					float dist1 = getManhattanDistance(meanLocs[0], tloc, nsnps);
					float dist2 = getManhattanDistance(meanLocs[1], tloc, nsnps);
					if (dist1 < dist2 && isInCluster1[trial][t] == false) {
						isInCluster1[trial][t] = true;
						noChanges = false;
					} else if (dist1 > dist2 && isInCluster1[trial][t] == true){
						isInCluster1[trial][t] = false;
						noChanges = false;
					}
				}

				if (noChanges) break;
			}
			
			if (badclusters == true) {
				System.out.println("Trial " + trial + ": bad clustering, no distance could be calculated");
			} else {
				float distanceBetweenClusters = getManhattanDistance(meanLocs[0], meanLocs[1], nsnps);
				if (distanceBetweenClusters > maxDistance) {
					maxDistance = distanceBetweenClusters;
					bestTrial = trial;
				}
				System.out.println("Trial " + trial + ": distance between clusters = " + distanceBetweenClusters);
			}
		}

		
		//make alignments based on the clusters
		boolean[] isInCluster2 = new boolean[ntaxa];
		for (int t = 0; t < ntaxa; t++) isInCluster2[t] = !isInCluster1[bestTrial][t];
		IdGroup id1 = IdGroupUtils.idGroupSubset(myAlignment.getIdGroup(), isInCluster1[bestTrial]);
		IdGroup id2 = IdGroupUtils.idGroupSubset(myAlignment.getIdGroup(), isInCluster2);
		
		Alignment a1 = FilterAlignment.getInstance(myAlignment, id1);
		Alignment a2 = FilterAlignment.getInstance(myAlignment, id2);
		
		return new Alignment[]{a1, a2};
	}
	
	public static float[] snpsAsFloatVector(BitSet[] alleles, int nsnps) {
		float[] result = new float[nsnps];
		for (int s = 0; s < nsnps; s++) {
			if (alleles[0].fastGet(s)) {
				result[s] = 2;
				if (alleles[1].fastGet(s)) result[s] = 1;
			} else {
				if (alleles[1].fastGet(s)) result[s] = 0;
//				else result[s] = -1;
				else result[s] = 1;
			}
		}
		return result;
	}
	
	public static float getManhattanDistance(float[] loc, float[] t, int nsnps) {
		float d = 0;
		int nsites = 0;
		for (int s = 0; s < nsnps; s++) {
			if (loc[s] >= 0 && t[s] >= 0) {
				d += Math.abs(loc[s] - t[s]);
				nsites++;
			}
		}
		return d / nsites;
	}
	
	public static float[] getMeanLocation(float[] loc, int[] size, float[] t, boolean add, int nsnps) {
		float[] result = new float[nsnps];
		if (add) {
			for (int s = 0; s < nsnps; s++) {
				if (t[s] >= 0) {
					if (size[s] > 0) {
						result[s] = (loc[s] * size[s] + t[s]) / ((float) (size[s] + 1));
						size[s]++;
					} else {
						result[s] = t[s];
						size[s] = 1;
					}
				} else {
					result[s] = loc[s];
				}
			}
		} else {
			for (int s = 0; s < nsnps; s++) {
				if (t[s] >= 0) {
					if (size[s] > 1) {
						result[s] = (loc[s] * size[s] - t[s]) / ((float) (size[s] - 1));
						size[s]--;
					} else if (size[s] == 1){
						result[s] = 0;
						size[s] = 0;
					}
				} else {
					result[s] = loc[s];
				}
			}
		}
		
		return result;
	}
	
	public static float[] getMeanLocation(float[][] locs) {
		int nsites = locs[0].length;
		int ntaxa = locs.length;
		float[] result = new float[nsites];
		
		for (int s = 0; s < nsites; s++) {
			int count = 0;
			float sum = 0;
			for (int t = 0; t < ntaxa; t++) {
				if (!Float.isNaN(locs[t][s])) {
					count++;
					sum += locs[t][s];
				}
				if (count > 0) result[s] = sum / count;
				else result[s] = Float.NaN;
			}
		}
		return result;
	}
	
	public static void printAlleleStats(Alignment a, String name) {
		int monoCount = 0;
		int polyCount = 0;
		int[] binCount = new int[21];
		int nsites = a.getSiteCount();
		for (int s = 0; s < nsites; s++) {
			if (a.getMajorAlleleFrequency(s) > 0.75) monoCount++;
			else {
				polyCount++;
				int bin = (int) Math.floor(20 * a.getMajorAlleleFrequency(s));
				binCount[bin]++;
			}
		}
		System.out.println(name);
		System.out.println("mono count = " + monoCount + ", poly count = " + polyCount);
		System.out.print("bins: ");
		for (int i = 0; i < 20; i++) System.out.print(" " + binCount[i]);
		System.out.println();
		System.out.println();
	}
	
	public static void mergeNonconsensusFiles(String dir, String match, String outfileName) {
		File[] mergeFiles = filterFiles(dir, match);
		
		int nfiles = mergeFiles.length;
		String[] colLabel = new String[nfiles];
		HashMap<String, String[]> taxonMap = new HashMap<String, String[]>();
		
		String input;
		String[] info;
		for (int f = 0; f < nfiles; f++) {
			System.out.println("processing" + mergeFiles[f].getName());
			try {
				BufferedReader br = new BufferedReader(new FileReader(mergeFiles[f]));
				br.readLine();
				input = br.readLine();
				info = tab.split(input);
				colLabel[f] = info[1];
				while (input != null) {
					info = tab.split(input);
					String[] values = taxonMap.get(info[0]);
					if (values == null) {
						values = new String[nfiles];
						for (int i = 0; i < nfiles; i++) values[i] = "";
						taxonMap.put(info[0], values);
					}
					values[f] = info[2];
					input =  br.readLine();
				}
				
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		File outfile = new File(dir, outfileName);
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
			StringBuilder sb = new StringBuilder("Taxon");
			for (int i = 0; i < nfiles; i++) {
				sb.append("\t").append(colLabel[i]);
			}
			bw.write(sb.toString());
			bw.newLine();
			LinkedList<String> taxaList = new LinkedList<String>(taxonMap.keySet());
			Collections.sort(taxaList);
			for (String taxon : taxaList) {
				sb = new StringBuilder(taxon);
				String[] values = taxonMap.get(taxon);
				for (int f = 0; f < nfiles; f++) sb.append("\t").append(values[f]);
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public static File[] filterFiles(String dir, String match) {
		File matchdir = new File(dir);
		final String pattern = new String(match);
		File[] filteredFiles = matchdir.listFiles(new FilenameFilter() {
			
			@Override
			public boolean accept(File dir, String name) {
				if (name.matches(pattern)) return true;
				return false;
			}
		});
		return filteredFiles;
	}
	
	public static void mergeFiles(File[] mergeFiles, int idcol, int datacol, int[] colOrder, String outfile) {
		int nfiles = mergeFiles.length;
		String input;
		String[] info;
		if (colOrder == null) colOrder = new int[]{0,2,3,4,5,6,7,8,9,1};
		HashMap<String, String[]> taxonMap = new HashMap<String, String[]>();
		
		for (int f = 0; f < nfiles; f++) {
			System.out.println("processing" + mergeFiles[f].getName());
			try {
				BufferedReader br = new BufferedReader(new FileReader(mergeFiles[f]));
				br.readLine();
				input = br.readLine();
				info = tab.split(input);
				while (input != null) {
					info = tab.split(input);
					String[] values = taxonMap.get(info[idcol]);
					if (values == null) {
						values = new String[nfiles];
						for (int i = 0; i < nfiles; i++) values[i] = "";
						taxonMap.put(info[idcol], values);
					}
					values[f] = info[datacol];
					input =  br.readLine();
				}
				
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
			StringBuilder sb = new StringBuilder("Taxon");
			for (int i = 0; i < nfiles; i++) {
				sb.append("\t").append(i + 1);
			}
			sb.append("\t").append("average");
			bw.write(sb.toString());
			bw.newLine();
			LinkedList<String> taxaList = new LinkedList<String>(taxonMap.keySet());
			Collections.sort(taxaList);
			for (String taxon : taxaList) {
				sb = new StringBuilder(taxon);
				String[] values = taxonMap.get(taxon);
				double sum = 0;
				double count = 0;
				for (int f = 0; f < nfiles; f++) {
					sb.append("\t").append(values[f]);
					try {
						sum += Double.parseDouble(values[f]);
						count++;
					} catch (NumberFormatException e) {
						//do nothing
					}
				}
				sb.append("\t").append(sum/count);
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public static void imputeLinkageMarkers(double interval, boolean hapmapFormat, String origsnpFile, String snpfilePattern, String outfilePattern) {
		String[] nuc = new String[]{"A","M","C"};
		
		HashMap<Byte, String> byteToNumberString = new HashMap<Byte, String>();
		byteToNumberString.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("A"), "0");
		byteToNumberString.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("M"), "1");
		byteToNumberString.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("C"), "2");
		
		HashMap<Byte, Double> byteToNumber = new HashMap<Byte, Double>();
		byteToNumber.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("A"), 0.0);
		byteToNumber.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("M"), 1.0);
		byteToNumber.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("C"), 2.0);
		
		Pattern tab = Pattern.compile("\t");
		BufferedWriter bw = null;
		byte missingByte = NucleotideAlignmentConstants.getNucleotideDiploidByte("NN");
		int chromosome = 1;
		String outFilename = String.format(outfilePattern, interval);
		//String snpFilename = "/Volumes/Macintosh HD 2/data/namgbs/genos_1217/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_chr" + chromosome + ".hmp.abhv2.impute5to3stateHMM.txt";
		String snpFilename = origsnpFile;
//		if (hapmapFormat) outFilename = "/Volumes/Macintosh HD 2/data/namgbs/genos_1217/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_hmp_" + interval + "cmsnps.hmp.txt";
//		else outFilename = "/Volumes/Macintosh HD 2/data/namgbs/genos_1217/NAM282_20111217_scv10mF8maf002_mgs_E1pLD5_hmp_" + interval + "cmsnps.txt";
		
		AGPMap agpmap = new AGPMap();
		
		int ntaxa = 0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(snpFilename));
			bw = new BufferedWriter(new FileWriter(outFilename));
			String header = br.readLine(); 
			String[] info = tab.split(header);
			int ncol = info.length;
			ntaxa = ncol - 11;
			
			if (hapmapFormat) bw.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
			else bw.write("Snp\tallele\tchr\tpos\tcm");
			for (int t = 11; t < ncol; t++) {
				bw.write("\t");
				bw.write(info[t]);
			}
			bw.newLine();
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		String[] family = new String[]{"Z001","Z002","Z003","Z004","Z005","Z006","Z007","Z008","Z009","Z010","Z011","Z012","Z013","Z014","Z015","Z016","Z018","Z019","Z020","Z021","Z022","Z023","Z024","Z025","Z026"};
		
		//impute data for each chromosome
		for (int chr = 1; chr <=10; chr++) {
			String chrstr = Integer.toString(chr);
			for (int fam = 0; fam < 25; fam++) {
				System.out.println("Imputing data for chromosome " + chr + ", family " + family[fam] + ".");
//				snpFilename = "/Volumes/Macintosh HD 2/results/recombination study/nam/final.Panzea/namibm.combined.hapmap.f.05r.5.chr" + chr + ".family."+ family[fam] + "parents.hmp.txt";
				snpFilename = String.format(snpfilePattern, chr, family[fam]);
				Alignment a = ImportUtils.readFromHapmap(snpFilename, true, null);
				int nsnps = a.getSiteCount();
				
				double startgenpos = agpmap.getCmFromPosition(chr, a.getPositionInLocus(0));
				//round up to nearest interval
				startgenpos = ((double) (Math.ceil(startgenpos / interval))) * interval;
				
				double endgenpos = agpmap.getCmFromPosition(chr, a.getPositionInLocus(nsnps - 1));
				//round down to nearest interval
				endgenpos = ((double)(Math.floor(endgenpos / interval))) * interval;

				
				int leftflank = 0;
				int rightflank = 0;
				try {
					for (double curpos = startgenpos; curpos <= endgenpos; curpos += interval) {
						int physpos = agpmap.getPositionFromCm(chr, curpos);
						String physposString = Integer.toString(physpos);
						String genpos = Double.toString(curpos);
						bw.write("S_");
						bw.write(physposString);
						bw.write("\timputed\t");
						bw.write(chrstr);
						bw.write("\t");
						bw.write(physposString);
						bw.write("\t");
						bw.write(genpos);
						if (hapmapFormat) bw.write("\tNA\tNA\tNA\tNA\tNA\tNA");
						while (physpos > a.getPositionInLocus(rightflank)) rightflank++; 
						leftflank = rightflank - 1;
						
						if (hapmapFormat) {
							for (int t = 0; t < ntaxa; t++) {
								bw.write("\t");
								int leftndx = leftflank;
								int rightndx = rightflank;
								while (a.getBase(t, leftndx) ==  missingByte && leftndx > 0) leftndx--;
								while (a.getBase(t, rightndx) ==  missingByte && rightndx < nsnps - 1) rightndx++;
								byte leftByte = a.getBase(t, leftndx);
								byte rightByte = a.getBase(t, rightndx);
								if (leftByte ==  missingByte) {
									if (rightByte == missingByte) bw.write("N");
									else bw.write(NucleotideAlignmentConstants.getNucleotideIUPAC(rightByte));
								}
								else if (rightByte ==  missingByte) bw.write(NucleotideAlignmentConstants.getNucleotideIUPAC(leftByte));
								else if (a.getBase(t, leftndx) == a.getBase(t, rightndx)) bw.write(NucleotideAlignmentConstants.getNucleotideIUPAC(leftByte));
								else bw.write("N"); 
							}
						} else {
							for (int t = 0; t < ntaxa; t++) {
								bw.write("\t");
								int leftndx = leftflank;
								int rightndx = rightflank;
								while (a.getBase(t, leftndx) ==  missingByte && leftndx > 0) leftndx--;
								while (a.getBase(t, rightndx) ==  missingByte && rightndx < nsnps - 1) rightndx++;
								byte leftByte = a.getBase(t, leftndx);
								byte rightByte = a.getBase(t, rightndx);
								if (leftByte ==  missingByte) {
									if (rightByte == missingByte) bw.write("-");
									else bw.write(byteToNumberString.get(NucleotideAlignmentConstants.getNucleotideIUPAC(rightByte)));
								}
								else if (rightByte ==  missingByte) bw.write(byteToNumberString.get(NucleotideAlignmentConstants.getNucleotideIUPAC(leftByte)));
								else if (a.getBase(t, leftndx) == a.getBase(t, rightndx)) bw.write(byteToNumberString.get(NucleotideAlignmentConstants.getNucleotideIUPAC(rightByte)));
								else {
									double leftval = byteToNumber.get(NucleotideAlignmentConstants.getNucleotideIUPAC(leftByte));
									double rightval = byteToNumber.get(NucleotideAlignmentConstants.getNucleotideIUPAC(rightByte));
									int leftpos = a.getPositionInLocus(leftndx);
									int rightpos = a.getPositionInLocus(rightndx);
									double pd = ((double) (physpos - leftpos)) / ((double) (rightpos - leftpos));
									double thisval = leftval * (1 - pd) + rightval * pd;
									bw.write(Double.toString(thisval));
								}; 
							}
						}
						bw.newLine();
					}
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(-1);
				}

			}
			

			
			
		}
		
		try {
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		
		System.out.println("Finished imputing markers.");
	}
	
	//use this for imputation of July final build results
	//modified for later builds
	public static void imputeLinkageMarkersAcrossFamilies(double interval, boolean hapmapFormat, boolean excludeTaxa) {
		class ImputedSnp {
			int physicalPos;
			double geneticPos;
			StringBuilder sb = new StringBuilder();
		}
		
		LinkedList<String> excludeList = getListOfTaxa("/Volumes/Macintosh HD 2/results/recombination study/nam/final.Panzea.consolidated.B/Nam.exclude.release.1.txt");
		String[] nuc = new String[]{"A","M","C"};

		Pattern tab = Pattern.compile("\t");
		byte missingByte = NucleotideAlignmentConstants.getNucleotideDiploidByte("NN");

		AGPMap agpmap = new AGPMap();

		//impute data for each chromosome
		for (int chr = 1; chr <=10; chr++) {

//			File snpfiledir = new File("/Volumes/Macintosh HD 2/results/recombination study/nam/final.Panzea.consolidated.B");
//			File snpfiledir = new File("/Volumes/Macintosh HD 2/results/recombination study/nam/final.Panzea.consolidated");
			File snpfiledir = new File("/Volumes/Macintosh HD 2/data/zea/build2.6/nam/imputed.var.plusfounders");
			final String chrname = "chr" + chr + ".family";
			File[] snpFiles = snpfiledir.listFiles(new FilenameFilter() {
				@Override
				public boolean accept(File dir, String name) {
//					if (name.startsWith("nam.consolidated") && name.contains(chrname)) return true;
					if (name.startsWith("USNAM_build2.6_imputed_var") && name.contains(chrname)) return true;
					return false;
				}
			});

			String chrstr = Integer.toString(chr);

			//Get the minimum and maximum physical positions for all families
			int minpos = Integer.MAX_VALUE;
			int maxpos = 0;
			for (File snpfile : snpFiles) {
				try {
					BufferedReader br = new BufferedReader(new FileReader(snpfile));
					br.readLine();
					String input = br.readLine();
					String[] info = tab.split(input, 5);
					int firstpos = Integer.parseInt(info[3]);
					String testInput;
					while ((testInput = br.readLine()) != null) {
						input = testInput;
					}
					info = tab.split(input, 5);
					int lastpos = Integer.parseInt(info[3]);

					minpos = Math.min(minpos, firstpos);
					maxpos = Math.max(maxpos, lastpos);
					br.close();
				} catch(IOException e) {
					e.printStackTrace();
					System.exit(-1);
				}
			}

			double startgenpos = agpmap.getCmFromPosition(chr, minpos);
			
			//round up to nearest interval
			startgenpos = ((double) (Math.ceil(startgenpos / interval))) * interval;
			int nIntervals = (int) ((agpmap.getCmFromPosition(chr, maxpos) - startgenpos) / interval);
			int nImputedSnps = nIntervals + 1;
			LinkedList<ImputedSnp> snpList = new LinkedList<ImputedSnp>();

			double curpos = startgenpos;
			for (int i = 0; i < nImputedSnps; i++) {
				ImputedSnp snp = new ImputedSnp();
				snp.geneticPos = curpos;
				curpos += interval;
				snp.physicalPos = agpmap.getPositionFromCm(chr, curpos);
				snpList.add(snp);
			}
			
			StringBuilder taxaHeader = new StringBuilder();
			for (File snpfile : snpFiles) {
				System.out.println("Imputing data for " + snpfile.getName() + ".");
				Alignment a = ImportUtils.readFromHapmap(snpfile.getPath(), true, null);
				
				boolean b73isA = isB73HaplotypeA(a);
				HashMap<Byte, String> byteToNumberString = new HashMap<Byte, String>();
				HashMap<Byte, Double> byteToNumber = new HashMap<Byte, Double>();
				HashMap<Byte, String> byteToNucleotide = new HashMap<Byte, String>();
				if (b73isA) {
					byteToNumberString.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("AA"), "0");
					byteToNumberString.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("AC"), "1");
					byteToNumberString.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("CA"), "1");
					byteToNumberString.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("CC"), "2");

					byteToNumber.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("AA"), 0.0);
					byteToNumber.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("AC"), 1.0);
					byteToNumber.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("CA"), 1.0);
					byteToNumber.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("CC"), 2.0);
					
					byteToNucleotide.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("AA"), "A");
					byteToNucleotide.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("AC"), "M");
					byteToNucleotide.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("CA"), "M");
					byteToNucleotide.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("CC"), "C");

				} else {
					byteToNumberString.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("AA"), "2");
					byteToNumberString.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("AC"), "1");
					byteToNumberString.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("CA"), "1");
					byteToNumberString.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("CC"), "0");

					byteToNumber.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("AA"), 2.0);
					byteToNumber.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("AC"), 1.0);
					byteToNumber.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("CA"), 1.0);
					byteToNumber.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("CC"), 0.0);
					
					byteToNucleotide.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("AA"), "C");
					byteToNucleotide.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("AC"), "M");
					byteToNucleotide.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("CA"), "M");
					byteToNucleotide.put(NucleotideAlignmentConstants.getNucleotideDiploidByte("CC"), "A");
				}
				
				int nsnps = a.getSiteCount();
				int ntaxa = a.getSequenceCount();
				int leftflank = 0;
				int rightflank = 0;
				
				for (int t = 0; t < ntaxa; t++) {
//					if (!a.getTaxaName(t).startsWith(excludeTaxon) && !excludeList.contains(a.getTaxaName(t))) taxaHeader.append("\t").append(a.getFullTaxaName(t));
//					if (a.getTaxaName(t).startsWith("Z0")) taxaHeader.append("\t").append(a.getFullTaxaName(t));
					if (useTaxon(a.getTaxaName(t), excludeList)) taxaHeader.append("\t").append(a.getFullTaxaName(t));
				}

				for (ImputedSnp isnp : snpList) {
					while (rightflank < nsnps && isnp.physicalPos > a.getPositionInLocus(rightflank)) rightflank++;
//					System.out.println("rightflank= " + rightflank + ", snp physicalPos= " + isnp.physicalPos + ", position of rightflank= " + a.getPositionInLocus(rightflank)); //debug
					leftflank = rightflank - 1;
//					System.out.println("leftflank= " + leftflank + ", snp physicalPos= " + isnp.physicalPos + ", position of leftflank= " + a.getPositionInLocus(leftflank)); //debug

					if (hapmapFormat) {
						for (int t = 0; t < ntaxa; t++) {
//							if (a.getTaxaName(t).startsWith(excludeTaxon)) continue;
//							if (excludeList.contains(a.getTaxaName(t))) continue;
//							if (!a.getTaxaName(t).startsWith("Z0")) continue;
							if (!useTaxon(a.getTaxaName(t), excludeList)) continue;
							isnp.sb.append("\t");
							byte leftByte, rightByte;
							
							if (leftflank < 0) leftByte = missingByte;
							else {
								int leftndx = leftflank;
								while (a.getBase(t, leftndx) ==  missingByte && leftndx > 0) leftndx--;
								leftByte = a.getBase(t, leftndx);
							}
							
							if (rightflank > nsnps - 1) rightByte = missingByte;
							else {
								int rightndx = rightflank;
								while (a.getBase(t, rightndx) ==  missingByte && rightndx < nsnps - 1) rightndx++;
								rightByte = a.getBase(t, rightndx);
							}
							if (leftByte ==  missingByte) {
								if (rightByte == missingByte) isnp.sb.append("N");
								else isnp.sb.append(byteToNucleotide.get(rightByte));
							} else if (rightByte ==  missingByte) isnp.sb.append(byteToNucleotide.get(leftByte));
							else if (leftByte == rightByte) isnp.sb.append(byteToNucleotide.get(leftByte));
							else isnp.sb.append("N"); 
						}
					} else {
						for (int t = 0; t < ntaxa; t++) {
//							if (a.getTaxaName(t).startsWith(excludeTaxon)) continue;
//							if (excludeList.contains(a.getTaxaName(t))) continue;
//							if (!a.getTaxaName(t).startsWith("Z0")) continue;
							if (!useTaxon(a.getTaxaName(t), excludeList)) continue;

							isnp.sb.append("\t");
							byte leftByte, rightByte;
							
							int leftndx = leftflank;
							if (leftflank < 0) leftByte = missingByte;
							else {
								while (a.getBase(t, leftndx) ==  missingByte && leftndx > 0) leftndx--;
								leftByte = a.getBase(t, leftndx);
							}
							
							int rightndx = rightflank;
							if (rightflank > nsnps - 1) rightByte = missingByte;
							else {
								while (a.getBase(t, rightndx) ==  missingByte && rightndx < nsnps - 1) rightndx++;
								rightByte = a.getBase(t, rightndx);
							}
							
							if (leftByte ==  missingByte) {
								if (rightByte == missingByte) isnp.sb.append("-");
								else isnp.sb.append(byteToNumberString.get(rightByte));
							}
							else if (rightByte ==  missingByte) isnp.sb.append(byteToNumberString.get(leftByte));
							else if (leftByte == rightByte) isnp.sb.append(byteToNumberString.get(rightByte));
							else {
								double leftval = byteToNumber.get(leftByte);
								double rightval = byteToNumber.get(rightByte);
								int leftpos = a.getPositionInLocus(leftndx);
								int rightpos = a.getPositionInLocus(rightndx);
								double pd = ((double) (isnp.physicalPos - leftpos)) / ((double) (rightpos - leftpos));
								double thisval = leftval * (1 - pd) + rightval * pd;
								isnp.sb.append(Double.toString(thisval));
							}
//							System.out.println("leftflank = " + leftflank + ", rightflank = " + rightflank + ", leftndx = " + leftndx + ", rightndx = " + rightndx + ", leftbyte = " + NucleotideAlignmentConstants.getNucleotideIUPAC(leftByte) + ", rightbyte = " + NucleotideAlignmentConstants.getNucleotideIUPAC(rightByte) + ", imputed value = " + isnp.sb.charAt(isnp.sb.length() - 1));
							
						}
					}
					
				}
			}
			
			//write out the chromosome file here
			File outfile;
//			if (hapmapFormat) outfile = new File(snpfiledir, "imputedMarkers.release.1.chr" + chrstr +"." + interval +"cm.final.Panzea.consolidated.B.hmp.txt");
//			else outfile = new File(snpfiledir, "imputedMarkers.release.1.chr" + chrstr +"." + interval +"cm.final.Panzea.consolidated.B.txt");
//			if (hapmapFormat) outfile = new File(snpfiledir, "imputedMarkers.release.1.chr" + chrstr +"." + interval +"cm.final.Panzea.consolidated.hmp.txt");
//			else outfile = new File(snpfiledir, "imputedMarkers.release.1.chr" + chrstr +"." + interval +"cm.final.Panzea.consolidated.txt");
//			if (hapmapFormat) outfile = new File(snpfiledir, "imputed.Markers/imputedMarkers.chr" + chrstr +"." + interval +"cm.final.Panzea.consolidated.B.hmp.txt");
//			else outfile = new File(snpfiledir, "imputed.Markers/imputedMarkers.chr" + chrstr +"." + interval +"cm.final.Panzea.consolidated.B.txt");
			if (hapmapFormat) outfile = new File(snpfiledir, "imputed.Markers/imputedMarkers.chr" + chrstr +"." + interval +"cm.USNAM2.6.imputed.var.hmp.txt");
			else outfile = new File(snpfiledir, "imputed.Markers/imputedMarkers.chr" + chrstr +"." + interval +"cm.USNAM2.6.imputed.var.txt");

			try{
				BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
				if (hapmapFormat) {
					bw.write("rs#\talleles\tchrom\tpos\tcm\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
				} else {
					bw.write("Snp\tallele\tchr\tpos\tcm");
				}
				bw.write(taxaHeader.toString());
				bw.write("\n");
				for (ImputedSnp isnp : snpList) {
					bw.write(String.format("S%d_%d\tNA\t", chr, isnp.physicalPos));
					bw.write(chrstr);
					bw.write("\t");
					bw.write(Integer.toString(isnp.physicalPos));
					bw.write("\t");
					bw.write(Double.toString(isnp.geneticPos));
					if (hapmapFormat) bw.write("\tNA\tNA\tNA\tNA\tNA\tNA");
					bw.write(isnp.sb.toString());
					bw.write("\n");
				}
				bw.close();
			} catch(IOException e) {

			}

		}
		
	}
	
	public static boolean useTaxon(String name, LinkedList<String> excludelist){
		if (excludelist != null) {
			if (name.startsWith("Z0") & !excludelist.contains(name)) return true;
		} else {
			if (name.startsWith("Z0")) return true;
		}
		
		return false;
	}
	
	public static boolean isB73HaplotypeA(Alignment a) {
		IdGroup ids = a.getIdGroup();
		int ndx = ids.whichIdNumber("B73(PI550473):MRG:2:250027110");
//		int ndx = ids.whichIdNumber("B73(PI550473)");
		int nsnps = a.getSiteCount();
		HashMap<Byte, Integer> alleleCounts = new HashMap<Byte, Integer>();
		for (int s = 0; s < nsnps; s++) {
			Byte allele = a.getBase(ndx, s);
			Integer count = alleleCounts.get(allele);
			if (count == null) alleleCounts.put(allele, 1);
			else alleleCounts.put(allele, 1 + count);
		}
		
//		System.out.println("B73 allele counts:");
		Byte bestAllele = -1;
		int maxCount = 0;
		for (Byte b:alleleCounts.keySet()) {
			int count = alleleCounts.get(b);
//			System.out.println(NucleotideAlignmentConstants.getNucleotideIUPAC(b) + ", " + count);
			if (count > maxCount) {
				maxCount = count;
				bestAllele = b;
			}
		}
		
		if (bestAllele.byteValue() == NucleotideAlignmentConstants.getNucleotideDiploidByte('A')) return true;
		
		return false;
	}
	
	public static void imputeLinkageMarkersFrom1106(double interval) {
		//the input data
		String snpFilename = "/Volumes/Macintosh HD 2/data/namgbs/ImputedMarkerGenotypes_flowering_traits_092909.txt";
		String outFilename = "/Volumes/Macintosh HD 2/results/namgbs/jointlinkage/array1106/imputedMarkers.1106.allchr.txt";
		ArrayList<float[]> genotypes = new ArrayList<float[]>();
		ArrayList<String> taxanames = new ArrayList<String>();
		int nmarkers = 1106;
		int ntaxa;

		AGPMap agpmap = new AGPMap();
		
		BufferedWriter bw = null;

		try {
			BufferedReader br = new BufferedReader(new FileReader(snpFilename));
			br.readLine();
			String input;
			while ((input = br.readLine()) != null) {
				String[] info = tab.split(input);
				float[] geno = new float[nmarkers];
				for (int i = 0; i < nmarkers; i++) geno[i] = Float.parseFloat(info[i+5]);
				genotypes.add(geno);
				taxanames.add(info[0]);
			}
			br.close();

			bw = new BufferedWriter(new FileWriter(outFilename));
			bw.write("Snp\tallele\tchr\tpos\tcm");
			for (String taxon:taxanames) {
				bw.write("\t");
				bw.write(taxon);
			}
			bw.write("\n");
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

		//impute data for each chromosome
		ntaxa = taxanames.size();
		for (int chr = 1; chr <=10; chr++) {
			String chrstr = Integer.toString(chr);
			double startgenpos = agpmap.getFirstGeneticPosition(chr);
			//round down to nearest interval
			startgenpos = ((double) (Math.floor(startgenpos / interval))) * interval;

			double endgenpos = agpmap.getLastGeneticPosition(chr);
			//round up to nearest interval
			endgenpos = ((double)(Math.ceil(endgenpos / interval))) * interval;

			try {
				for (double curpos = startgenpos; curpos <= endgenpos; curpos += interval) {
					int physpos = agpmap.getPositionFromCm(chr, curpos);
					String physposString = Integer.toString(physpos);
					String genpos = Double.toString(curpos);
					bw.write("S");
					bw.write(chrstr);
					bw.write("_");
					bw.write(physposString);
					bw.write("\timputed\t");
					bw.write(chrstr);
					bw.write("\t");
					bw.write(physposString);
					bw.write("\t");
					bw.write(genpos);
					int[] flanks = agpmap.getFlankingMarkerIndices(chr, curpos);
					
					for (int t = 0; t < ntaxa; t++) {
						bw.write("\t");
						double val;
						
						if (flanks[0] < 0) val = genotypes.get(t)[flanks[1]];
						else if(flanks[1] >= nmarkers) val = genotypes.get(t)[flanks[0]];
						else if (flanks[0] == flanks[1]) val = genotypes.get(t)[flanks[0]]; 
						else {
							double pd = (curpos - agpmap.getGeneticPos(flanks[0])) / (agpmap.getGeneticPos(flanks[1]) - agpmap.getGeneticPos(flanks[0]));
							val = genotypes.get(t)[flanks[0]] * (1 - pd) + genotypes.get(t)[flanks[1]] * pd;
						}
						bw.write(Double.toString(val));
						
					}
					
					bw.newLine();
				}
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}

		}

		try {
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
		}

		System.out.println("Finished imputing markers.");


	}
	
	public static LinkedList<String> getListOfTaxa(String filename) {
		LinkedList<String> taxaList = new LinkedList<String>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String taxon;
			while((taxon = br.readLine()) != null) taxaList.add(taxon);
			br.close();
		} catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		
		return taxaList;
	}
	
	
}

