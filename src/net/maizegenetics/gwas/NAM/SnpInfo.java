package net.maizegenetics.gwas.NAM;

public class SnpInfo {
	final int chromosome;
	final int pos;
	final String allele;
	final double[] genotype;
	final double F;
	final double p;
	
	SnpInfo(int chromosome, int pos, String allele, double[] genotype, double F, double p) {
		this.chromosome = chromosome;
		this.pos = pos;
		this.allele = allele;
		this.genotype = genotype;
		this.F = F;
		this.p = p;
	}

}
