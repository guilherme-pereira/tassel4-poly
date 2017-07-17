package net.maizegenetics.gwas.NAM;

public class SnpDataB73Ref extends SnpData {

	public SnpDataB73Ref(int chromosome) {
		super(chromosome);
	}
	
	public SnpDataB73Ref(FileNames files) {
		super(files);
	}

	@Override
	public double[] getGenotype() {
		int count = 0;
		double[] geno = new double[26];
		double b73 = Double.parseDouble(parsedLine[11]);
		for (int i = 12; i < 38; i++) {
			geno[count++] = Math.abs(b73 - Double.parseDouble(parsedLine[i]));
		}
		
		return geno;
	}

	@Override
	public boolean[] getBooleanGenotype() {
		int count = 0;
		boolean[] geno = new boolean[26];
		String b73 = parsedLine[11];
		for (int i = 12; i < 38; i++) {
			geno[count++] = (parsedLine[i].equals(b73));
		}
		return geno;
	}

}
