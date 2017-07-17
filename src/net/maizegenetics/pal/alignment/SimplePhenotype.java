package net.maizegenetics.pal.alignment;

import java.util.List;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;

public class SimplePhenotype extends AbstractPhenotype {
	private DoubleMatrix2D data;
	
	public SimplePhenotype(IdGroup taxa, List<Trait> traits, DoubleMatrix2D data) {
		super(taxa, traits);
		this.data = data;
	}

	public SimplePhenotype(IdGroup taxa, List<Trait> traits, double[][] data) {
		super(taxa, traits);
		this.data = DoubleFactory2D.dense.make(data);
	}
	
	public SimplePhenotype(IdGroup taxa, List<Trait> traits) {
		super(taxa, traits);
		this.data = DoubleFactory2D.dense.make(getNumberOfTaxa(), getNumberOfTraits(), Double.NaN);
	}
	
	public double getData(int taxon, int trait) {
		return data.get(taxon, trait);
	}

	public double getData(Identifier taxon, Trait trait) {
		return getData(whichTaxon(taxon), whichTrait(trait));
	}

	public double[][] getData() {
		return data.toArray();
	}

	public void setData(int taxon, int trait, double value) {
		data.set(taxon, trait, value);
	}

	public void setData(Identifier taxon, Trait trait, double value) {
		data.set(whichTaxon(taxon), whichTrait(trait), value);
	}

}
