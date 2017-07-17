package net.maizegenetics.gwas.NAM;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;


public class SnpData {
	int chromosome;
	BufferedReader br = null;
	String[] parsedLine;
	int numberOfSnps;
	FileNames files = null;
	int leftmarker;
	double relativeDistance;
	
	public final static String[] popnames = new String[]{"B97","CML103","CML228","CML247","CML277","CML322","CML333","CML52R","CML69","HP301","IL14H",
			"KI11","KI3","KY21","M162W","M37W","MO17","MO18W","MS71","NC350","NC358","OH43","OH7B","P39","TX303","TZI8"};
		

	public SnpData(int chromosome) {
		this.chromosome = chromosome;
	}
	
	public SnpData(FileNames files) {
		this.chromosome = files.chromosome;
		this.files = files;
		findTotalSnpNumber();
		reset();
	}
	
	public void findTotalSnpNumber() {
		reset();
		numberOfSnps = 0;
		try {
			while (br.readLine() != null) numberOfSnps++;
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
	}
	
	public boolean next() {
	
		if (br == null) return false;
		try {
			String inLine = br.readLine();
			if (inLine == null) {
				br.close();
//				dis.close();
				return false;
			}
			else {
				parsedLine = inLine.split("\t");
//				leftmarker = dis.readInt();
//				relativeDistance = (double) dis.readFloat();
			}
			
		} catch (IOException e) {
			e.printStackTrace();
			return false;
		}
		
//		if (IsB73Alt())  return next();
		
		return true;
	}
	
	public double[] getGenotype() {
		int count = 0;
		double[] geno = new double[26];
		for (int i = 12; i < 38; i++) {
			geno[count++] = Double.parseDouble(parsedLine[i]);
		}
		return geno;
	}
	
	//equals true for the B73 allele
	public boolean[] getBooleanGenotype() {
		int count = 0;
		boolean[] geno = new boolean[26];
		for (int i = 12; i < 38; i++) {
			geno[count++] = (parsedLine[i].equals("0"));
		}
		return geno;
	}
	
	public int getPosition() {
		int pos = Integer.parseInt(parsedLine[3]);
		return Integer.parseInt(parsedLine[3]);
	}
	
	public void reset() {

		try {
			if (br != null) br.close();
			br = new BufferedReader(new FileReader(files.snps));
			br.readLine();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public boolean IsB73Alt() {
		char alt = parsedLine[1].charAt(2);
		if (parsedLine[11].charAt(0) == alt) return true;
		return false;
	}
	
	public int getNumberOfSnps() {return numberOfSnps;}
	
	public String getAllele() {return parsedLine[1];}
	
	public SnpData getCopy() {
		SnpData newData = new SnpData(chromosome);
		newData.files = files;
		newData.numberOfSnps = numberOfSnps;
		return newData;
	}
	
}
