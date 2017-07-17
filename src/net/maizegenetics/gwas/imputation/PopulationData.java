package net.maizegenetics.gwas.imputation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.regex.Pattern;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.BitSet;

public class PopulationData {
	public String name;
	public ArrayList<String> members;
	public BitSet snpIndex;
	public String parent1;
	public String parent2;
	public double contribution1;
	public double contribution2;
	public int Fgen;
	public double inbredCoef = -1;
	public Alignment original;
	public Alignment imputed;
	public byte[] alleleA;
	public byte[] alleleC;
	public boolean useSubpopulations;
	public ArrayList<IdGroup> subpopulationGroups;
	public ArrayList<BitSet> subpoplationSiteIndex;

	/**
	 * @param popFilename	the name of the file containing pedigree information for a group of populations
	 * @return	a HashMap of family names (keys) and associated PopulationData objects (values) containing information for the families in the pedigree file
	 */
	public static ArrayList<PopulationData> readPedigreeFile(String popFilename) {
		Pattern tab = Pattern.compile("\t");
		LinkedHashMap<String, PopulationData> familyMap = new LinkedHashMap<String, PopulationData>();
		String input = "";
		try {
			BufferedReader br = new BufferedReader(new FileReader(popFilename));
			br.readLine();
			while ((input = br.readLine()) != null) {
				String[] info = tab.split(input);
				if (info[0].length() > 0 && !info[0].equalsIgnoreCase("NA")) {
					PopulationData family = familyMap.get(info[0]);
					if (family == null) {
						family = new PopulationData ();
						family.name = info[0];
						family.members = new ArrayList<String>();
						family.members.add(info[1]);
						family.parent1 = info[2];
						family.parent2 = info[3];
						family.contribution1 = Double.parseDouble(info[4]);
						family.contribution2 = Double.parseDouble(info[5]);
						try {
							family.inbredCoef = Double.parseDouble(info[6]);
						} catch (Exception e) {}
						
						familyMap.put(info[0], family);
					}
					else family.members.add(info[1]);
				}
			}
			br.close();
		} catch (IOException e) {
			System.out.println("at: " + input);
			e.printStackTrace();
			System.exit(-1);
		} catch (NumberFormatException nfe) {
			System.out.println("at: " + input);
			nfe.printStackTrace();
			System.exit(-1);
		}
		return new ArrayList<PopulationData>(familyMap.values());
	}
}