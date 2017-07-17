package net.maizegenetics.pal.alignment;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeSet;

import net.maizegenetics.jGLiM.BasicLevel;

import net.maizegenetics.pal.ids.Identifier;

public class MarkerPhenotypeAdapter {

    protected MarkerPhenotype markerpheno;
    protected Phenotype pheno;
    protected Alignment align;
    protected int numberOfPhenotypes;
    protected int numberOfFactorTraits;
    protected int numberOfCovariates;
    protected int numberOfMarkers;
    protected int numberOfMarkersFromAlignment;
    protected int numberOfMarkersFromPhenotype;
    protected ArrayList<String> phenotypeNames;
    protected ArrayList<String> covariateNames;
    protected ArrayList<String> factorNames;
    protected ArrayList<ArrayList<PhenotypeInfo>> blockList;
    protected HashMap<String, HashMap<BasicLevel, Integer>> covariateTraitMap;
    protected HashMap<String, HashMap<BasicLevel, Integer>> factorTraitMap;
    protected ArrayList<Integer> markerIndex;

    public MarkerPhenotypeAdapter(Phenotype aPhenotype) {
        markerpheno = null;
        align = null;
        pheno = aPhenotype;
        init();
    }

    public MarkerPhenotypeAdapter(MarkerPhenotype aMarkerPhenotype) {
        markerpheno = aMarkerPhenotype;
        align = aMarkerPhenotype.getAlignment();
        pheno = aMarkerPhenotype.getPhenotype();
        init();
    }

    protected void init() {
        //get info from alignment
        if (align == null) {
            numberOfMarkersFromAlignment = 0;
        } else {
            numberOfMarkersFromAlignment = align.getSiteCount();
        }

        //now deal with the phenotype
        //are Trait Factors consistent?
        if (!areTraitFactorsConsistent()) {
            throw new IllegalArgumentException("Some traits have inconsistent headers. No analysis will be performed");
        }

        //the columns will be one for each factor, one for each unique trait name
        //get the unique trait names
        TreeSet<String> datanameSet = new TreeSet<String>();
        covariateTraitMap = new HashMap<String, HashMap<BasicLevel, Integer>>();
        factorTraitMap = new HashMap<String, HashMap<BasicLevel, Integer>>();
        markerIndex = new ArrayList<Integer>();
        int nTraits = pheno.getNumberOfTraits();

        for (int i = 0; i < nTraits; i++) {
            Trait trait = pheno.getTrait(i);
            if (trait.getType().equals(Trait.TYPE_DATA)) {
                datanameSet.add(trait.getName());
            }
            if (trait.getType().equals(Trait.TYPE_COVARIATE)) {
                String name = trait.getName();
                HashMap<BasicLevel, Integer> covmap = covariateTraitMap.get(name);
                if (covmap == null) {
                    covmap = new HashMap<BasicLevel, Integer>();
                    covariateTraitMap.put(name, covmap);
                }
                covmap.put(getTraitLevel(trait), i);
            }
            if (trait.getType().equals(Trait.TYPE_FACTOR)) {
                String name = trait.getName();
                HashMap<BasicLevel, Integer> facmap = factorTraitMap.get(name);
                if (facmap == null) {
                    facmap = new HashMap<BasicLevel, Integer>();
                    factorTraitMap.put(name, facmap);
                }
                facmap.put(getTraitLevel(trait), i);
            }
            if (trait.getType().equals(Trait.TYPE_MARKER)) {
                markerIndex.add(i);
            }
        }

        numberOfPhenotypes = datanameSet.size();
        phenotypeNames = new ArrayList<String>(datanameSet);
        numberOfCovariates = covariateTraitMap.size();
        covariateNames = new ArrayList<String>(covariateTraitMap.keySet());
        Collections.sort(covariateNames);
        numberOfFactorTraits = factorTraitMap.size();
        factorNames = new ArrayList<String>(factorTraitMap.keySet());
        Collections.sort(factorNames);
        numberOfMarkersFromPhenotype = markerIndex.size();
        numberOfMarkers = numberOfMarkersFromAlignment + numberOfMarkersFromPhenotype;

        //create a list of blocks for each phenotype
        blockList = new ArrayList<ArrayList<PhenotypeInfo>>();
        for (int p = 0; p < numberOfPhenotypes; p++) {
            blockList.add(new ArrayList<PhenotypeInfo>());
        }

        //Each block is a combination of the phenotype factors
        for (int i = 0; i < nTraits; i++) {
            Trait trait = pheno.getTrait(i);
            if (trait.getType().equals(Trait.TYPE_DATA)) {
                int ndx = Collections.binarySearch(phenotypeNames, trait.getName());
                BasicLevel level = getTraitLevel(trait);
                blockList.get(ndx).add(new PhenotypeInfo(i, level, trait.getName()));
            }
        }

        String msg = checkFactorCovariateEnvironments();
        if (msg.length() > 0) {
            throw new IllegalArgumentException("No analysis will be run:\n" + msg);
        }

    }

    public boolean areTraitFactorsConsistent() {
        //traits should all have the same factors or no factors
        if (pheno.getNumberOfFactors() == 0) {
            return true;
        }
        String[] factors = pheno.getFactorNameCopy();
        for (Trait trait : pheno.getTraits()) {
            ArrayList<String> factorlist = trait.getFactorNames();
            if (factorlist != null && factorlist.size() > 0) {
                if (factorlist.size() != factors.length) {
                    return false;
                }
                for (int i = 0; i < factors.length; i++) {
                    if (!factorlist.contains(factors[i])) {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    public String checkFactorCovariateEnvironments() {
        //factor columns should either have no phenotype factors or the ones they have should match at least one data column
        //start by creating a list of factor levels that exist for phenotypes
        StringBuilder msg = new StringBuilder();
        HashSet<BasicLevel> levelSet = new HashSet<BasicLevel>();
        for (ArrayList<PhenotypeInfo> phenoList : blockList) {
            for (PhenotypeInfo pi : phenoList) {
                levelSet.add(pi.level);
            }
        }

        //check factors
        for (Map.Entry<String, HashMap<BasicLevel, Integer>> entry : factorTraitMap.entrySet()) {
            boolean consistent = true;
            for (Map.Entry<BasicLevel, Integer> ent : entry.getValue().entrySet()) {
                BasicLevel level = ent.getKey();
                if (level != null && (level.getNumberOfSublevels() > 0 && !levelSet.contains(level))) {
                    consistent = false;
                }
            }
            if (!consistent) {
                String factorname = entry.getKey();
                msg.append("The factor " + factorname + " has environments that are inconsistent with the data.\n");
            }
        }

        //check covariates
        for (Map.Entry<String, HashMap<BasicLevel, Integer>> entry : covariateTraitMap.entrySet()) {
            boolean consistent = true;
            for (Map.Entry<BasicLevel, Integer> ent : entry.getValue().entrySet()) {
                BasicLevel level = ent.getKey();
                if (level != null && (level.getNumberOfSublevels() > 0 && !levelSet.contains(level))) {
                    consistent = false;
                }
            }
            if (!consistent) {
                String covname = entry.getKey();
                msg.append("The covariate " + covname + " has environments that are inconsistent with the data.\n");
            }
        }

        return msg.toString();
    }

    public BasicLevel getTraitLevel(Trait trait) {
        String[] factorName = pheno.getFactorNameCopy();
        int nFactors = factorName.length;
        if (trait.getNumberOfFactors() != nFactors) {
            return null;
        }
        String[] level = new String[nFactors];
        for (int i = 0; i < nFactors; i++) {
            level[i] = trait.getFactorValue(factorName[i]);
        }
        return new BasicLevel(level);
    }

    public int getNumberOfRows(int phenotype) {
        int blocks = blockList.get(phenotype).size();
        return pheno.getNumberOfTaxa() * blocks;
    }

    public int getNumberOfBlocks(int phenotype) {
        return blockList.get(phenotype).size();
    }

    public int getNumberOfPhenotypes() {
        return numberOfPhenotypes;
    }

    public double[] getPhenotypeValues(int i) {
        double[] phenotype = new double[getNumberOfRows(i)];
        int ntaxa = pheno.getNumberOfTaxa();
        int block = 0;
        for (PhenotypeInfo info : blockList.get(i)) {
            for (int t = 0; t < ntaxa; t++) {
                phenotype[block * ntaxa + t] = pheno.getData(t, info.index);
            }
            block++;
        }
        return phenotype;
    }

    public boolean[] getMissingPhenotypes(int phenotype) {
        return MarkerPhenotypeAdapterUtils.whichAreMissing(getPhenotypeValues(phenotype));
    }

    public String getPhenotypeName(int i) {
        return phenotypeNames.get(i);
    }

    public int getNumberOfFactors() {
        return numberOfFactorTraits + pheno.getNumberOfFactors();
    }

    public String[] getFactorValues(int phenotype, int factor) {
        int nrows = getNumberOfRows(phenotype);
        String[] factorvalue = new String[nrows];
        int ntaxa = pheno.getNumberOfTaxa();
        int numberOfPhenoFactors = pheno.getNumberOfFactors();

        if (factor < numberOfPhenoFactors) {
            factorvalue = new String[nrows];
            int block = 0;
            for (PhenotypeInfo info : blockList.get(phenotype)) {
                for (int t = 0; t < ntaxa; t++) {
                    factorvalue[block * ntaxa + t] = (String) info.level.getSublevel(factor);
                }
                block++;
            }
        } else {
            factor -= numberOfPhenoFactors;
            factorvalue = new String[nrows];
            int block = 0;
            for (PhenotypeInfo info : blockList.get(phenotype)) {
                Integer ndx = factorTraitMap.get(factorNames.get(factor)).get(null);
                if (ndx == null) {
                    ndx = factorTraitMap.get(factorNames.get(factor)).get(info.level);
                }
                String[] factorlevel = pheno.getTrait(ndx).getLevelLabels();
                for (int t = 0; t < ntaxa; t++) {
                    factorvalue[block * ntaxa + t] = factorlevel[(int) pheno.getData(t, ndx)];
                }
                block++;
            }
        }

        return factorvalue;
    }

    public boolean[] getMissingFactors(int phenotype, int factor) {
        return MarkerPhenotypeAdapterUtils.whichAreMissing(getFactorValues(phenotype, factor));
    }

    public String getFactorName(int i) {
        if (i < pheno.getNumberOfFactors()) {
            return pheno.getFactorName(i);
        } else {
            return factorNames.get(i - pheno.getNumberOfFactors());
        }
    }

    public int getNumberOfCovariates() {
        return numberOfCovariates;
    }

    public double[] getCovariateValues(int phenotype, int covariate) {
        double[] cov = new double[getNumberOfRows(phenotype)];
        int ntaxa = pheno.getNumberOfTaxa();
        int block = 0;
        for (PhenotypeInfo info : blockList.get(phenotype)) {
            Integer ndx = covariateTraitMap.get(covariateNames.get(covariate)).get(null);
            if (ndx == null) {
                ndx = covariateTraitMap.get(covariateNames.get(covariate)).get(info.level);
            }
            for (int t = 0; t < ntaxa; t++) {
                if (ndx != null) {
                    cov[block * ntaxa + t] = pheno.getData(t, ndx);
                } else {
                    cov[block * ntaxa + t] = Double.NaN;
                }
            }
            block++;
        }
        return cov;
    }

    public boolean[] getMissingCovariates(int phenotype, int covariate) {
        return MarkerPhenotypeAdapterUtils.whichAreMissing(getCovariateValues(phenotype, covariate));
    }

    public String getCovariateName(int i) {
        return covariateNames.get(i);
    }

    public int getNumberOfMarkers() {
        return numberOfMarkers;
    }

    public boolean isMarkerDiscrete(int i) {
        if (i < numberOfMarkersFromAlignment) {
            return true;
        }
        i -= numberOfMarkersFromAlignment;
        return pheno.getTrait(markerIndex.get(i)).isDiscrete();
    }

    public String getMarkerName(int i) {
        if (i < numberOfMarkersFromAlignment) {
            return align.getSNPID(i);
        }
        int t = markerIndex.get(i - numberOfMarkersFromAlignment);
        return pheno.getTrait(t).getName();
    }

    public String getMarkerChromosome(int marker) {
        if (marker < numberOfMarkersFromAlignment) {
            return align.getLocus(marker).getChromosomeName();
        }
        int t = markerIndex.get(marker - numberOfMarkersFromAlignment);
        Object prop = pheno.getTrait(t).getProperty(Trait.PROP_CHROMOSOME);
        if (prop == null) {
            return "";
        }
        return prop.toString();
    }

    public double getMarkerChromosomePosition(int marker) {
        if (marker < numberOfMarkersFromAlignment) {
            return Double.NaN;
        }
        int t = markerIndex.get(marker - numberOfMarkersFromAlignment);
        Object chrpos = pheno.getTrait(t).getProperty(Trait.PROP_POSITION);
        if (chrpos == null) {
            return Double.NaN;
        }
        return ((Double) chrpos).doubleValue();
    }

    public String getLocusName(int marker) {
        if (marker < numberOfMarkersFromAlignment) {
            return align.getLocusName(marker);
        }
        int t = markerIndex.get(marker - numberOfMarkersFromAlignment);
        Object locusName = pheno.getTrait(t).getProperty(Trait.PROP_LOCUS);
        if (locusName == null) {
            return "Unknown";
        }
        return locusName.toString();
    }

    public int getLocusPosition(int marker) {
        if (marker < numberOfMarkersFromAlignment) {
            return align.getPositionInLocus(marker);
        }
        int t = markerIndex.get(marker - numberOfMarkersFromAlignment);
        Object prop = pheno.getTrait(t).getProperty(Trait.PROP_LOCUS_POSITION);
        if (prop == null || !(prop instanceof Integer)) {
            return marker;
        }
        return ((Integer) prop).intValue();
    }

    public Object[] getMarkerValue(int phenotype, int marker) {
        int nrows = getNumberOfRows(phenotype);
        int nblocks = getNumberOfBlocks(phenotype);
        int ntaxa = pheno.getNumberOfTaxa();
        if (marker < numberOfMarkersFromAlignment) {
            String[] values = new String[nrows];
            for (int b = 0; b < nblocks; b++) {
                for (int t = 0; t < ntaxa; t++) {
                    values[b * ntaxa + t] = markerpheno.getAlignment().getBaseAsString(t, marker);
                }
            }
            return values;
        } else {
            int ndx = markerIndex.get(marker - numberOfMarkersFromAlignment);
            if (isMarkerDiscrete(marker - numberOfMarkersFromAlignment)) {
                String[] values = new String[nrows];
                for (int b = 0; b < nblocks; b++) {
                    String[] labels = pheno.getTrait(ndx).getLevelLabels();
                    for (int t = 0; t < ntaxa; t++) {
                        double val = pheno.getData(t, ndx);
                        if (Double.isNaN(val)) {
                            values[b * ntaxa + t] = Phenotype.MISSING_STRING;
                        } else {
                            values[b * ntaxa + t] = labels[(int) val];
                        }
                    }
                }
                return values;
            } else {
                Double[] values = new Double[nrows];
                for (int b = 0; b < nblocks; b++) {
                    for (int t = 0; t < ntaxa; t++) {
                        values[b * ntaxa + t] = pheno.getData(t, markerIndex.get(marker));
                    }
                }
                return values;
            }
        }
    }

    public boolean[] getMissingMarkers(int phenotype, int marker) {
        Object[] values = getMarkerValue(phenotype, marker);
        return MarkerPhenotypeAdapterUtils.whichAreMissing(values);
    }

    public Identifier[] getTaxa(int phenotype) {
        int nrows = getNumberOfRows(phenotype);
        Identifier[] names = new Identifier[nrows];
        int ntaxa = pheno.getNumberOfTaxa();
        int nblocks = getNumberOfBlocks(phenotype);
        for (int b = 0; b < nblocks; b++) {
            for (int t = 0; t < ntaxa; t++) {
                names[b * ntaxa + t] = pheno.getTaxon(t);
            }
        }
        return names;
    }

    public boolean hasMarkerNames() {
        if (numberOfMarkersFromPhenotype > 0) {
            return true;
        }
        if (numberOfMarkersFromAlignment > 0) {
            if (align.getSNPID(0) != null) {
                return true;
            }
        }

        return false;
    }

    class PhenotypeInfo {

        int index;
        BasicLevel level;
        String name;

        PhenotypeInfo(int index, BasicLevel level, String name) {
            this.index = index;
            this.level = level;
            this.name = name;
        }
    }
}
