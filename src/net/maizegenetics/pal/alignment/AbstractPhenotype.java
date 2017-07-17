package net.maizegenetics.pal.alignment;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.LinkedHashSet;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.pal.ids.Identifier;

public abstract class AbstractPhenotype implements Phenotype {
	private String[] factors;
	private boolean[] useFactor;
	private IdGroup taxa;
	private List<Trait> traitList;

	public AbstractPhenotype(IdGroup taxa, List<Trait> traitList) {
		this.taxa = taxa;
		this.traitList = traitList;
		init();
	}
	
	public AbstractPhenotype(IdGroup taxa, List<Trait> traitList, String[] factors, boolean[] useFactor) {
		this.taxa = taxa;
		this.traitList = traitList;
		this.factors = factors;
		this.useFactor = useFactor;
	}
	
	private void init() {
		LinkedHashSet<String> factorList = new LinkedHashSet<String>();
		for (Trait trait : traitList) {
			ArrayList<String> theseFactors = trait.getFactorNames();
			if (theseFactors != null) factorList.addAll(theseFactors);
		}
		
		int n = factorList.size();
		factors = new String[n];
		useFactor = new boolean[n];
		if (n > 0) factorList.toArray(factors);
		for (int i = 0; i < n; i++) useFactor[i] = true;
	}
	
	public String getActiveFactorName(int factor) {
		int count = 0;
		int active = -1;
		int n = factors.length;
		
		while (count < n) {
			if (useFactor[count]) active++;
			if (factor == active) return factors[count];
			count++;
		}
	
		return null;
	}

	public String[] getActiveFactorNames() {
		int n = getNumberOfActiveFactors();
		String[] names = new String[n];
		int count = 0;
		int activeCount = 0;
		for (String name : factors) {
			if (useFactor[count++]) names[activeCount++] = name;
		}
		return names;
	}

	public String getFactorName(int factor) {
		return factors[factor];
	}

	public String[] getFactorNameCopy() {
		String[] names = new String[factors.length];
		System.arraycopy(factors, 0, names, 0, names.length);
		return names;
	}

	public int getNumberOfActiveFactors() {
		int count = 0;
		for (boolean active : useFactor) if (active) count++;
		return count;
	}

	public int getNumberOfFactors() {
		return factors.length;
	}

	public int getNumberOfTaxa() {
		return taxa.getIdCount();
	}

	public int getNumberOfTraits() {
		return traitList.size();
	}

	public IdGroup getTaxa() {
		return taxa;
	}

	public Identifier getTaxon(int taxon) {
		return taxa.getIdentifier(taxon);
	}

	public Trait getTrait(int trait) {
		return traitList.get(trait);
	}

	public List<Trait> getTraits() {
		return traitList;
	}

	public boolean isFactorActive(int factor) {
		return useFactor[factor];
	}

	public void setActiveFactor(int factor, boolean active) {
		useFactor[factor] = active;
	}

	public int whichTaxon(Identifier taxon) {
		return taxa.whichIdNumber(taxon);
	}

	public int whichTrait(Trait trait) {
		return traitList.indexOf(trait);
	}
	
	public static List<Trait> copyTraitsFromPhenotype(Phenotype phenotype) {
		List<Trait> newTraitList = new ArrayList<Trait>();
		for (Trait t : phenotype.getTraits()) newTraitList.add(t);
		return newTraitList;
	}

	@Override
	public void report(PrintWriter out) {
		//out.print("Data Table: \n");
		//out.println(getRowCount() + " taxa \n" + getColumnCount() + " columns (including taxa)");
	}

	@Override
	public int getColumnCount() {
		return traitList.size() + 1;
	}

	@Override
	public int getElementCount() {
		return getColumnCount() * getRowCount();
	}

	@Override
	public int getRowCount() {
		return taxa.getIdCount();
	}

	@Override
	public Object[] getTableColumnNames() {
		String[] names = new String[getColumnCount()];
		int nTraits = getNumberOfTraits();
		names[0] = "Taxa";
		for (int i = 0; i < nTraits; i++) {
			names[i + 1] = traitList.get(i).toString();
		}
		return names;
	}

	@Override
	public String getTableTitle() {
		String title = "Phenotypes";
		return title;
	}
	@Override
	public Object[] getRow(int row) {
		int n = getColumnCount();
		Object[] rowvalues = new Object[n];
		rowvalues[0] = getTaxon(row);
		for (int i = 1; i < n; i++) {
			Trait trait = getTrait(i - 1);
			if (trait.isDiscrete()) {
				String[] labels = (String[]) trait.getProperty(Trait.PROP_LEVEL_LABELS);
				double dblval = getData(row, i - 1);
				if (Double.isNaN(dblval)) rowvalues[i] = "?";
				else rowvalues[i] = labels[(int) dblval];
			}
			else rowvalues[i] = getData(row, i - 1);
		}
		return rowvalues;
	}

	@Override
	public Object[][] getTableData() {
		int n = getRowCount();
		Object[][] tabledata = new Object[n][];
		for (int i = 0; i < n; i++) tabledata[i] = getRow(i);
		return tabledata;
	}

	@Override
	public Object[][] getTableData(int start, int end) {
		int n = end - start + 1;
		Object[][] tabledata = new Object[n][];
		for (int i = start; i <= end; i++) tabledata[i] = getRow(i);
		return tabledata;
	}

	@Override
	public String[] getRowNames() {
		return IdGroupUtils.getNames(getTaxa());
	}


	@Override
	public Object getValueAt(int row, int col) {
		if (col == 0) return getTaxon(row);
		col --;
		Trait trait = traitList.get(col);
		if (trait.isDiscrete()) {
			double dblval = getData(row, col);
			if (Double.isNaN(dblval)) return "?";
			else return trait.getLevelLabel((int) dblval);
		}
		else {
			return new Double(getData(row, col));
		}
	}

	
}
