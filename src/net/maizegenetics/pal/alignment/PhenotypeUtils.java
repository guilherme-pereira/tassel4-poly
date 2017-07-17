package net.maizegenetics.pal.alignment;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

public class PhenotypeUtils {

	//prevents instantiation of this utility class
	private PhenotypeUtils() {
		
	}
	
	/**
	 * @param pheno	a Phenotype
	 * @param traitType	a trait type 
	 * @return	a List of all the traits in the Phenotype with a trait type of traitType. 
	 */
	public static List<Trait> getListofType(Phenotype pheno, String traitType) {
		ArrayList<Trait> traitList = new ArrayList<Trait>();
		int ntraits = pheno.getNumberOfTraits();
		for (int t = 0; t < ntraits; t++) {
			Trait trait = pheno.getTrait(t);
			if (trait.getType().equals(traitType)) traitList.add(trait);
		}
		return traitList;
	}
	
	/**
	 * @param pheno a Phenotype
	 * @return true if the Phenotype markers are all discrete or all continuous, false otherwise. 
	 * If the Phenotype has no markers, returns true.
	 */
	public static boolean areMarkersConsistentForDiscreteness(Phenotype pheno) {
		int ntraits = pheno.getNumberOfTraits();
		ArrayList<Boolean> discreteList = new ArrayList<Boolean>();
		for (int t = 0; t < ntraits; t++) {
			Trait trait = pheno.getTrait(t);
			if (trait.getType().equals(Trait.TYPE_MARKER)) discreteList.add(trait.isDiscrete());
		}
		
		int nmarkers = discreteList.size();
		Iterator<Boolean> it = discreteList.iterator();
		Boolean discrete = it.next();
		boolean consistent = true;
		while (it.hasNext() && consistent) {
			if (!discrete.equals(it.next())) consistent = false;
		}
		
		return consistent;
	}
	
	/**
	 * @param pheno a Phenotype
	 * @return true if all of the Phenotypes marker alleles have the same ploidy, false otherwise.
	 * If the Phenotype has no markers, returns true.
	 */
	public static boolean areMarkersConsistentForPloidy(Phenotype pheno) {
		Pattern colon = Pattern.compile(":");
		List<Trait> markerList = getListofType(pheno, Trait.TYPE_MARKER);
		if (!markerList.get(0).isDiscrete()) return true;
		String firstlabel = markerList.get(0).getLevelLabel(0);
		int nalleles = colon.split(firstlabel).length;
		
		boolean consistent = true;
		Iterator<Trait> it = markerList.iterator();
		while (it.hasNext() && consistent) {
			for (String label : it.next().getLevelLabels()) {
				if (nalleles != colon.split(label).length) consistent = false;
			}
		}
		return consistent;
	}
	
	/**
	 * used to save a Phenotype in a format that can be read by TASSEL
	 * @param pheno a Phenotype
	 * @param pw the PrintWriter to which the output will be written
	 */
	public static void saveAs(Phenotype pheno, PrintWriter pw) {
		StringBuilder sb;
		String sep = "\t";
		int nheaders = pheno.getNumberOfFactors();
		int ntraits = pheno.getNumberOfTraits();
		int ntaxa = pheno.getNumberOfTaxa();
		
		//output headers
		for (int h = 0; h < nheaders; h++) {
			sb = new StringBuilder("<Header name=");
			String factorName = pheno.getFactorName(h);
			sb.append(factorName).append(">");
			for (int t = 0; t < ntraits; t++) {
				sb.append(sep).append(pheno.getTrait(t).getFactorValue(factorName));
			}
			pw.println(sb.toString());
		}
		
		//output use
		sb = new StringBuilder("<Use>");
		for (int t = 0; t < ntraits; t++) {
			sb.append(sep).append(pheno.getTrait(t).getType());
		}
		pw.println(sb.toString());
		
		//output format
		sb = new StringBuilder("<Format>");
		for (int t = 0; t < ntraits; t++) {
			sb.append(sep);
			if (pheno.getTrait(t).isDiscrete()) sb.append("char");
			else sb.append("num");
		}
		pw.println(sb.toString());
		
		
		//output traitnames
		sb = new StringBuilder("<Trait>");
		for (int t = 0; t < ntraits; t++) {
			sb.append(sep).append(pheno.getTrait(t).getName());
		}
		pw.println(sb.toString());
		
		//output data
		int nrows = pheno.getRowCount();
		int ncol = pheno.getColumnCount();
		for (int r = 0; r < nrows; r++) {
			sb = new StringBuilder(pheno.getValueAt(r, 0).toString());
			for (int c = 1; c < ncol; c++) sb.append(sep).append(pheno.getValueAt(r, c));
			pw.println(sb.toString());
		}
		
	}
}
