package net.maizegenetics.gwas.modelfitter;

import java.util.ArrayList;

import net.maizegenetics.pal.alignment.Locus;

public class SNP {
	public String name;
	public Locus locus;
	public int position;
	ArrayList<Object> alleles;
	int index;
	
	public SNP(String name, Locus locus, int positionInLocus, int index) {
		this.name = name;
		this.locus = locus;
		position = positionInLocus;
		this.index = index;
	}
	
	public SNP() {
		
	}

	@Override
	public String toString() {
		return name;
	}
	
	
}
