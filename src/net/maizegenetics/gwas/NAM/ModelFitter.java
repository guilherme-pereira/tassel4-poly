package net.maizegenetics.gwas.NAM;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;


public class ModelFitter extends Thread {
	float[][] genotypes;
	int firstMarker;
	static AGPMap theAGPMap = null;
	double[][] residuals;
	String[] residualSamples;
	HashMap<String, Integer> sampleNameMap;
	int traitnumber;
	String thisTraitName;
	int totalSnps;
	int snpcount;
	FileNames files;
	
	SnpData snpdata;
	int maxmarker;
	ArrayList<String> pops;
	int[] popIndex;

	
	public static final int[] chromosomeLength = new int[]{300239041,
		234752839,
		230558137,
		247095508,
		216915529,
		169254300,
		170974187,
		174515299,
		152350485,
		149686046};

	public ModelFitter(int chromosome, FileNames files) {
		getAGPMap(files.agpmap);
		this.files = files;
	}
	
	public ModelFitter(FileNames parameters) {
		files = parameters;
		getAGPMap(files.agpmap);
	}
	
	public void run() {
		init();
		long start = System.currentTimeMillis();
		testSnps();
		long elapsed = System.currentTimeMillis() - start;
		System.out.println("Elapsed time for namgwas = " + elapsed);
	}
	
	public void init() {
		importResiduals(files.residuals);
		importNamMarkersforMap(files.namMarkersByChr);
		
		if (files.useB73asReference) snpdata = new SnpDataB73Ref(files);
		else if (files.randomizeSnpOrder) snpdata = new SnpDataRandomOrder(files);
		else snpdata = new SnpData(files);
		maxmarker = genotypes[0].length - 1;
		int nSamples = residualSamples.length;
		popIndex = new int[nSamples];
		pops = new ArrayList<String>();
		int count = 0;
		for (String sample : residualSamples) {
			int ipop = getPopulation(sample) - 1;
			pops.add(SnpData.popnames[ipop]);
			popIndex[count++] = ipop;
		}

	}
	
	public synchronized void getAGPMap(File mapfile) {
		if (theAGPMap == null) {
			theAGPMap = new AGPMap(mapfile);
		}
	}
	
	//reads the combined NAM + IBM data file for a chromosome
	//residuals must have been imported first so that the genotypes are properly matched to them
	protected void importNamMarkersforMap(File file) {
		Pattern tab = Pattern.compile("\t");
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			//the first line contains headers
			String input = br.readLine();
			String[] data = tab.split(input);
			firstMarker = Integer.parseInt(data[2].substring(1));

			int nMarkers = data.length - 1;
			int nTaxa = sampleNameMap.size();
			genotypes = new float[nTaxa][nMarkers];
			
			while ((input = br.readLine()) != null) {
				data = tab.split(input);
				Integer sampleNumber = sampleNameMap.get(data[0]);
				if (sampleNumber != null) {
					int n = sampleNumber.intValue();
					for (int j = 0; j < nMarkers; j++) {
						genotypes[n][j] = Float.parseFloat(data[j + 1]);
					}
				}
			}
			br.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	protected void importResiduals(File residualfile) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(residualfile));
			int nRows = 0;
			br.readLine();
			while (br.readLine() != null) nRows++;
			br.close();
			
			residuals = new double[nRows][10];
			residualSamples = new String[nRows];
			sampleNameMap = new HashMap<String, Integer>();
			br = new BufferedReader(new FileReader(residualfile));
			br.readLine();
			for (int i = 0; i < nRows; i++) {
				String[] parsedLine = br.readLine().split("\t");
				residualSamples[i] = parsedLine[0];
				sampleNameMap.put(residualSamples[i], new Integer(i));
				for (int c = 0; c < 10; c++) {
					residuals[i][c] = Double.parseDouble(parsedLine[c + 1]);
				}
			}
			br.close();
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

	}
	
	public void createProjectionData() {
		System.out.println("project data for chromosome " + files.chromosome);
		Pattern tab = Pattern.compile("\t");
		snpdata.reset();
		try {
			
			FileOutputStream fos = new FileOutputStream(files.projectedFile);
			BufferedOutputStream bos = new BufferedOutputStream(fos, 50000);
			DataOutputStream dos = new DataOutputStream(bos);
			
			int snpcount = 0;
			while (snpdata.next()) {
				if (snpcount % 1000 == 0) System.out.println("snp count = " + snpcount);
				snpcount++;
				int pos = snpdata.getPosition();
				boolean[] parents = snpdata.getBooleanGenotype();
				Object[][] markers = theAGPMap.getInterval(files.chromosome, pos);
				int left = AGPMap.getMarkerPosition(markers[0]);
				if (left == -1) left = 0;
				int right = AGPMap.getMarkerPosition(markers[1]);
				if (right == -1) right = chromosomeLength[files.chromosome - 1];
				
				//proportion of distance of snp between left and right markers
				float pd = ((float)(pos - left)) / ((float)(right - left));
				
				int leftmarker = AGPMap.getMarkerNumber(markers[0]);
				if (leftmarker == -1) leftmarker = 0; //telomere
				else leftmarker = leftmarker - firstMarker + 1;
				int rightmarker = AGPMap.getMarkerNumber(markers[1]);
				if (rightmarker == -1) rightmarker = maxmarker; //telomere
				else rightmarker = rightmarker - firstMarker + 1;

				int nSamples = residualSamples.length;
				float snpvalue;
				for (int i = 0; i < nSamples; i++) {
					if (parents[popIndex[i]]) {
						snpvalue = 0;
					} else {
						float leftval = (float) genotypes[i][leftmarker];
						float rightval = (float) genotypes[i][rightmarker];
						if (leftval == rightval) {
							snpvalue = leftval;
						} else {
							snpvalue = leftval * (1-pd) + rightval * pd;
						}
						
					}
					dos.writeFloat(snpvalue);
				}
			}
			
			dos.close();
			bos.close();
			fos.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	protected void testSnps() {
		//does nothing, needs to implemented in a subclass
		
	}
	
	public double[] projectSnp(boolean[] parents, int pos, int[] popIndex) {
		Object[][] markers = theAGPMap.getInterval(files.chromosome, pos);
		int left = AGPMap.getMarkerPosition(markers[0]);
		if (left == -1) left = 0;
		int right = AGPMap.getMarkerPosition(markers[1]);
		if (right == -1) right = chromosomeLength[files.chromosome - 1];
		
		//proportion of distance of snp between left and right markers
		double pd = ((double)(pos - left)) / ((double)(right - left));
		
		int leftmarker = AGPMap.getMarkerNumber(markers[0]);
		if (leftmarker == -1) leftmarker = 0; //telomere
		else leftmarker = leftmarker - firstMarker + 1;
		int rightmarker = AGPMap.getMarkerNumber(markers[1]);
		if (rightmarker == -1) rightmarker = maxmarker; //telomere
		else rightmarker = rightmarker - firstMarker + 1;

		int nSamples = residualSamples.length;
		double[] snpvalues = new double[nSamples];
		
		for (int i = 0; i < nSamples; i++) {
			if (parents[popIndex[i]]) {
				snpvalues[i] = 0;
			} else {
				double leftval = genotypes[i][leftmarker];
				double rightval = genotypes[i][rightmarker];
				if (leftval == rightval) {
					snpvalues[i] = leftval;
				} else {
					snpvalues[i] = leftval * (1-pd) + rightval * pd;
				}
				
			}
		}
		return snpvalues;

	}

	public double[] projectSnp(double[] parents, int pos, int[] popIndex) {
		int chromosome = files.chromosome;
		Object[][] markers = theAGPMap.getInterval(chromosome, pos);
		int left = AGPMap.getMarkerPosition(markers[0]);
		if (left == -1) left = 0;
		int right = AGPMap.getMarkerPosition(markers[1]);
		if (right == -1) right = chromosomeLength[chromosome - 1];
		
		//proportion of distance of snp between left and right markers
		double pd = ((double)(pos - left)) / ((double)(right - left));
		
		int leftmarker = AGPMap.getMarkerNumber(markers[0]);
		if (leftmarker == -1) leftmarker = 0; //telomere
		else leftmarker = leftmarker - firstMarker + 1;
		int rightmarker = AGPMap.getMarkerNumber(markers[1]);
		if (rightmarker == -1) rightmarker = maxmarker; //telomere
		else rightmarker = rightmarker - firstMarker + 1;

		int nSamples = residualSamples.length;
		double[] snpvalues = new double[nSamples];
		
		for (int i = 0; i < nSamples; i++) {
			double parent = parents[popIndex[i]];
			if (parent == 0) {
				snpvalues[i] = 0;
			} else {
				double leftval = genotypes[i][leftmarker];
				double rightval = genotypes[i][rightmarker];
				if (leftval == rightval) {
					snpvalues[i] = leftval * parent;
				} else {
					snpvalues[i] = (leftval * (1-pd) + rightval * pd) * parent;
				}
				
			}
		}
		return snpvalues;

	}

	public synchronized double[] projectSnp(boolean[] parents, int pos) {
		int chromosome = files.chromosome;
		Object[][] markers = theAGPMap.getInterval(chromosome, pos);
		int left = AGPMap.getMarkerPosition(markers[0]);
		if (left == -1) left = 0;
		int right = AGPMap.getMarkerPosition(markers[1]);
		if (right == -1) right = chromosomeLength[chromosome - 1];
		
		//proportion of distance of snp between left and right markers
		double pd = ((double)(pos - left)) / ((double)(right - left));
		
		int leftmarker = AGPMap.getMarkerNumber(markers[0]);
		if (leftmarker == -1) leftmarker = 0; //telomere
		else leftmarker = leftmarker - firstMarker + 1;
		int rightmarker = AGPMap.getMarkerNumber(markers[1]);
		if (rightmarker == -1) rightmarker = maxmarker; //telomere
		else rightmarker = rightmarker - firstMarker + 1;

		int nSamples = residualSamples.length;
		double[] snpvalues = new double[nSamples];
		
		for (int i = 0; i < nSamples; i++) {
			if (parents[popIndex[i]]) {
				snpvalues[i] = 0;
			} else {
				double leftval = genotypes[i][leftmarker];
				double rightval = genotypes[i][rightmarker];
				if (leftval == rightval) {
					snpvalues[i] = leftval;
				} else {
					snpvalues[i] = leftval * (1-pd) + rightval * pd;
				}
				
			}
		}
		return snpvalues;

	}

	public synchronized double[] projectSnp(double[] parents, int pos) {
		int chromosome = files.chromosome;
		Object[][] markers = theAGPMap.getInterval(chromosome, pos);
		int left = AGPMap.getMarkerPosition(markers[0]);
		if (left == -1) left = 0;
		int right = AGPMap.getMarkerPosition(markers[1]);
		if (right == -1) right = chromosomeLength[chromosome - 1];
		
		//proportion of distance of snp between left and right markers
		double pd = ((double)(pos - left)) / ((double)(right - left));
		
		int leftmarker = AGPMap.getMarkerNumber(markers[0]);
		if (leftmarker == -1) leftmarker = 0; //telomere
		else leftmarker = leftmarker - firstMarker + 1;
		int rightmarker = AGPMap.getMarkerNumber(markers[1]);
		if (rightmarker == -1) rightmarker = maxmarker; //telomere
		else rightmarker = rightmarker - firstMarker + 1;

		int nSamples = residualSamples.length;
		double[] snpvalues = new double[nSamples];
		
		for (int i = 0; i < nSamples; i++) {
			double parent = parents[popIndex[i]];
			if (parent == 0) {
				snpvalues[i] = 0;
			} else  {
				double leftval = genotypes[i][leftmarker];
				double rightval = genotypes[i][rightmarker];
				if (leftval == rightval) {
					snpvalues[i] = leftval * parent;
				} else {
					snpvalues[i] = (leftval * (1-pd) + rightval * pd) * parent;
				}
				
			}
		}
		return snpvalues;

	}

	protected BufferedWriter openOutputFile(String filename) {
		try {
			return new BufferedWriter(new FileWriter(filename));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

	protected BufferedWriter openOutputFile(File file) {
		try {
			return new BufferedWriter(new FileWriter(file));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}
	
	protected BufferedWriter openOutputFile(File file, boolean append) {
		try {
			return new BufferedWriter(new FileWriter(file, append));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}
	
	protected void writeToOutput(String out, BufferedWriter bw) {
		try {
			bw.write(out);
			bw.newLine();
			bw.flush();
		} catch (IOException e) {
			e.printStackTrace();
			try{bw.close();} catch (Exception e2) {}
			System.exit(-1);
		}
	}

	protected void closeOutput(BufferedWriter bw) {
		try {
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static int getPopulation(String sample) {
		if (sample.charAt(0) == 'M') return 17;
		return Integer.parseInt(sample.substring(1, 4));
	}
	
	public double getProgress() {
		return ((double) snpcount) / ((double) totalSnps);
	}

	public String getThisTraitName() {
		return thisTraitName;
	}

	public void setThisTraitName(String thisTraitName) {
		this.thisTraitName = thisTraitName;
	}
	
}
