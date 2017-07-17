package net.maizegenetics.gwas.imputation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.regex.Pattern;

public class TransitionProbabilityWithVariableRecombination extends TransitionProbability{
//	double averageRecombinationRate;
	int[] ratePosition = null;
	double[] rateAtPosition = null;
	int numberOfPositions = 0;
	int[][] transitionCounts;
	int numberOfTaxa;
	
	public TransitionProbabilityWithVariableRecombination(String chromosome) {
		readRateFile(chromosome);
		setInitialTransitionCounts(Integer.parseInt(chromosome));
	}

	@Override
	public void setNode(int node) {
		if (node <= 0) return;
		double rrr = getRelativeRecombinationRate(node);
		int n = transitionCounts.length;
		adjustedProbability = new double[n][n];
		int segmentLength = positions[node] - positions[node - 1];
		double m;
		for (int row = 0; row < n; row++) {
			double offdiagsum = 0;
			for (int col = 0; col < n; col++) {
				if (col != row) {
					adjustedProbability[row][col] = segmentLength * rrr * transitionCounts[row][col] / numberOfTaxa;
					offdiagsum += adjustedProbability[row][col];
				}
			}
			adjustedProbability[row][row] = 1 - offdiagsum;
		}
	}
	
//	private double getRelativeRecombinationRate(int node) {
//		int start = Arrays.binarySearch(ratePosition, positions[node - 1]);
//		int end = Arrays.binarySearch(ratePosition, positions[node]);
//		if (start < 0) start = -start - 2;
//		if (start < 0) start = 0;
//		if (end < 0) end = -end - 1;
//		if (end >= numberOfPositions) end = numberOfPositions - 1;
//		double sum = 0;
//		for (int i = start; i<=end; i++) sum += rateAtPosition[i];
//		
//		return sum/(end - start + 1);
//	}
	
	@Override
	public void setTransitionCounts(int[][] transitionCounts,
			int chromosomeLength, int numberOfTaxa) {
		
		super.setTransitionCounts(transitionCounts, chromosomeLength, numberOfTaxa);
		this.transitionCounts = transitionCounts;
		this.numberOfTaxa = numberOfTaxa;
	}

	private double getRelativeRecombinationRate(int node) {
		int index1 = Arrays.binarySearch(ratePosition, positions[node]);
		if (index1 < 0) {
			index1 = -index1 - 1;
		}
		if (index1 >= ratePosition.length) index1 = ratePosition.length - 1;
		int index0 = Arrays.binarySearch(ratePosition, positions[node - 1]);
		if (index0 < 0) {
			index0 = -index0 - 2;
		}
		if (index0 < 0 ) index0 = 0;
		
		double count = 0;
		double sum = 0;
		for (int i = index0; i <= index1; i++) {
			sum += Math.max(rateAtPosition[i], 0);
			count ++;
		}
		return sum/count;
	}
	
	private boolean readRateFile(String chromosome) {
		String filename = String.format("/Volumes/Macintosh HD 2/data/zea/build2.6/nam/imputed/xo.rate/xo.rate.chr%s.spar.5.txt", chromosome);
		Pattern tab = Pattern.compile("\t");
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			int nlines = 0;
			while (br.readLine() != null) nlines++;
			br.close();
			int npts = nlines; 
			ratePosition = new int[npts];
			rateAtPosition = new double[npts];
			br = new BufferedReader(new FileReader(filename));
			for (int i = 0; i < npts; i++) {
				String[] data = tab.split(br.readLine());
				ratePosition[i] = (int) Double.parseDouble(data[0]);
				rateAtPosition[i] = Double.parseDouble(data[1]);
			}
			br.close();
			numberOfPositions = npts;
			return true;
		} catch(IOException e) {
			e.printStackTrace();
			return false;
		}
	}
	
	public void setInitialTransitionCounts(int chr) {
		///these counts are from imputed US NAM population Z001, which works as a starting point
		numberOfTaxa = 195;
		switch(chr) {
		case 1:
			transitionCounts = new int[][]{
					{296522,32,82,28,291},
					{38,11049,18,10,34},
					{84,25,26051,29,69},
					{31,11,21,14385,43},
					{293,33,85,40,313513}};
			break;
		case 2:
			transitionCounts = new int[][]{
					{255745,41,90,22,231},
					{33,8755,25,6,29},
					{59,28,20773,28,63},
					{31,2,20,10326,25},
					{256,23,48,22,235366}				
			};
			break;
		case 3:
			transitionCounts = new int[][]{
					{218458,32,75,14,216},
					{28,9898,33,0,36},
					{62,36,15865,20,74},
					{18,3,21,5489,13},
					{207,29,62,19,165602}					
			};
			break;
		case 4:
			transitionCounts = new int[][]{
					{198184,17,35,36,217},
					{31,8092,9,7,21},
					{27,9,11230,11,40},
					{33,3,14,7087,20},
					{211,38,26,17,170986}					
			};
			break;
		case 5:
			transitionCounts = new int[][]{
					{249231,28,56,28,221},
					{33,8072,5,3,22},
					{50,7,15048,19,46},
					{20,3,14,10172,32},
					{210,24,45,20,195230}					
			};
			break;
		case 6:
			transitionCounts = new int[][]{
					{164187,13,40,19,185},
					{15,5598,7,0,15},
					{42,6,10974,13,52},
					{12,4,9,3325,18},
					{174,14,52,9,142522}					
			};
			break;
		case 7:
			transitionCounts = new int[][]{
					{200209,17,51,32,192},
					{27,8646,15,2,21},
					{52,12,11372,7,47},
					{26,5,6,6266,19},
					{189,25,46,15,152558}					
			};
			break;
		case 8:
			transitionCounts = new int[][]{
					{145044,13,62,12,191},
					{18,4657,9,3,9},
					{55,12,13521,17,53},
					{17,4,16,5892,18},
					{183,14,48,19,181770}					
			};
			break;
		case 9:
			transitionCounts = new int[][]{
					{141264,11,43,19,199},
					{11,3779,8,0,10},
					{37,5,11831,22,40},
					{18,2,10,7248,30},
					{198,12,44,22,166966}					
			};
			break;
		case 10:
			transitionCounts = new int[][]{
					{131018,28,30,16,135},
					{24,7285,14,3,16},
					{30,7,9196,2,44},
					{19,1,9,4372,19},
					{140,24,25,25,133253}					
			};
			break;
			
		}
	}
}
