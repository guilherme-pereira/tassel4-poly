package net.maizegenetics.baseplugins;

import java.awt.Frame;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import org.apache.log4j.Logger;

import net.maizegenetics.pal.alignment.FilterPhenotype;
import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

public class FilterTraitsPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(FilterTraitsPlugin.class);
	ArrayList<int[]> includeList;
	ArrayList<String[]> traitTypesList;
	
	public FilterTraitsPlugin(Frame parentFrame, boolean isInteractive) {
		super(parentFrame, isInteractive);
		
	}
	
	@Override
	public String getButtonName() {
		return "Traits";
	}

	@Override
	public ImageIcon getIcon() {
        URL imageURL = FilterTraitsPlugin.class.getResource("images/Filter.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
	}

	@Override
	public String getToolTipText() {
		return "Exclude traits or change properties";
	}

	@Override
	public DataSet performFunction(DataSet input) {
		List<Datum> data = input.getDataOfType(Phenotype.class);
		ArrayList<Datum> outputList = new ArrayList<Datum>();
		
		if (isInteractive()) {
            if (data.size() == 0) {
                JOptionPane.showMessageDialog(getParentFrame(), "No Phenotype data selected.");
            }
            includeList = new ArrayList<int[]>();
            traitTypesList = new ArrayList<String[]>();
			for (Datum datum : data) {
				FilterTraitsDialog ftd = new FilterTraitsDialog(getParentFrame(), (Phenotype) datum.getData());
                ftd.setLocationRelativeTo(getParentFrame());
				ftd.setVisible(true);
				includeList.add(ftd.getIncludedTraits());
				traitTypesList.add(ftd.getTraitTypes());
				ftd.dispose();
			}
		}

		int n = includeList.size();
		for (int i = 0; i < n; i++) {
			Datum datum = data.get(i);
			Phenotype pheno = (Phenotype) datum.getData();
			int[] included = includeList.get(i);
			if (included != null) {
                // This excludeLast flag is set if included has one
                // value of -1.  This to remove the last column of the
                // Phenotype data.  Mainly used by pipeline to
                // remove last column of population sturcture data.
                boolean excludeLast = false;
                if ((included.length == 1) && (included[0] == -1)) {
                    excludeLast = true;
                    included = new int[pheno.getNumberOfTraits() - 1];
                    for (int f = 0; f < (pheno.getNumberOfTraits() - 1); f++) {
                        included[f] = f;
                    }
                }
				FilterPhenotype filteredPhenotype = FilterPhenotype.getInstance(pheno, null, included);
				int ntraits = included.length;
                String[] types = null;
                if (excludeLast) {
                    types = new String[ntraits];
                    for (int f = 0; f < ntraits; f++) {
                        types[f] = pheno.getTrait(f).getType();
                    }
                } else {
				    types = traitTypesList.get(i);
                }
                
				for (int t = 0; t < ntraits; t++) {
					filteredPhenotype.getTrait(t).setType(types[included[t]]);
				}
				String name = "Filtered_" + datum.getName();
				StringWriter sw = new StringWriter();
				filteredPhenotype.report(new PrintWriter(sw));
				outputList.add(new Datum(name, filteredPhenotype, sw.toString()));
			}
		}
		
		
		DataSet ds = new DataSet(outputList, this);
		fireDataSetReturned(ds);
		return ds;
	}
	
	public void addIncludedTraits(int[] traitsToInclude) {
		includeList.add(traitsToInclude);
	}
	
	public void addTraitTypes(String[] types) {
		traitTypesList.add(types);
	}

	public void setIncludeList(ArrayList<int[]> includeList) {
		this.includeList = includeList;
	}

	public void setTraitTypesList(ArrayList<String[]> traitTypesList) {
		this.traitTypesList = traitTypesList;
	}
	
	
}


