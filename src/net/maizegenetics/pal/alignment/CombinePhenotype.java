package net.maizegenetics.pal.alignment;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.Map.Entry;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;

public class CombinePhenotype extends AbstractPhenotype {

    Phenotype[] phenotypes;
    int numberOfPhenotypes;
    boolean isUnion;
    int[][] phenotypeColumn; //first dimension is columns for this Phenotype, second dimension is phenotype index then offset
    int[][] phenotypeRow;	//first dimension is rows for this Phenotype, second dimension are rows in sub-phenotypes
    int totalColumns;
    int totalRows;

    private CombinePhenotype(Phenotype[] phenotypes, boolean isUnion, IdGroup taxaGroup, List<Trait> traitList, int[][] rowmap, int[][] columnmap) {
        super(taxaGroup, traitList);
        this.phenotypes = phenotypes;
        this.isUnion = isUnion;
        phenotypeRow = rowmap;
        phenotypeColumn = columnmap;
        totalColumns = phenotypeColumn.length;
        totalRows = phenotypeRow.length;
    }

    public static CombinePhenotype getInstance(Phenotype phenotype1, Phenotype phenotype2, boolean isUnion) {
        return getInstance(new Phenotype[]{phenotype1, phenotype2}, isUnion);
    }

    public static CombinePhenotype getInstance(Phenotype[] phenotypes, boolean isUnion) {
        Object[] rowinfo = mapRows(phenotypes, isUnion);
        Object[] colinfo = mapColumns(phenotypes, isUnion);

        int[][] rowmap = (int[][]) rowinfo[0];
        IdGroup taxa = (IdGroup) rowinfo[1];
        int[][] colmap = (int[][]) colinfo[0];
        ArrayList<Trait> traitList = (ArrayList<Trait>) colinfo[1];

        return new CombinePhenotype(phenotypes, isUnion, taxa, traitList, rowmap, colmap);
    }

    private static Object[] mapRows(Phenotype[] phenotypes, boolean isUnion) {
        TreeMap<Identifier, int[]> rowTreeMap = new TreeMap<Identifier, int[]>();
        int npheno = phenotypes.length;
        int pcount = 0;
        for (Phenotype pheno : phenotypes) {
            int n = pheno.getNumberOfTaxa();
            for (int t = 0; t < n; t++) {
                Identifier id = pheno.getTaxon(t);
                int[] rows = rowTreeMap.get(id);
                if (rows == null) {
                    rows = new int[npheno];
                    for (int i = 0; i < npheno; i++) {
                        rows[i] = -1;
                    }
                    rowTreeMap.put(id, rows);
                }
                rows[pcount] = t;
            }
            pcount++;
        }

        //if intersect, delete taxa that are not in all phenotypes
        Set<Entry<Identifier, int[]>> rowSet = rowTreeMap.entrySet();
        if (!isUnion) {
            Iterator<Entry<Identifier, int[]>> rit = rowSet.iterator();
            while (rit.hasNext()) {
                boolean notComplete = false;
                int[] rows = rit.next().getValue();
                for (int i : rows) {
                    if (i == -1) {
                        notComplete = true;
                        break;
                    }
                }
                if (notComplete) {
                    rit.remove();
                }
            }
        }

        int nrows = rowSet.size();
        int[][] rowmap = new int[nrows][];
        Identifier[] ids = new Identifier[nrows];
        int count = 0;
        for (Entry<Identifier, int[]> entry : rowSet) {
            rowmap[count] = entry.getValue();
            ids[count] = entry.getKey();
            count++;
        }

        return new Object[]{rowmap, new SimpleIdGroup(ids)};
    }

    private static Object[] mapColumns(Phenotype[] phenotypes, boolean isUnion) {
        int ncol = 0;
        for (Phenotype pheno : phenotypes) {
            ncol += pheno.getNumberOfTraits();
        }
        ArrayList<Trait> traitList = new ArrayList<Trait>(ncol);
        int[][] colmap = new int[ncol][2];
        int count1 = 0;
        int count2 = 0;
        for (Phenotype pheno : phenotypes) {
            int n = pheno.getNumberOfTraits();
            for (int i = 0; i < n; i++) {
                colmap[count2][0] = count1;
                colmap[count2][1] = i;
                count2++;
                traitList.add(Trait.getInstance(pheno.getTrait(i)));
            }
            count1++;
        }

        return new Object[]{colmap, traitList};
    }

    //implement Phenotype methods
    public double[][] getData() {
        int nrows = getNumberOfTraits();
        int ncols = getNumberOfTraits();
        double[][] result = new double[nrows][ncols];
        for (int r = 0; r < nrows; r++) {
            for (int c = 0; c < ncols; c++) {
                result[r][c] = getData(r, c);
            }
        }
        return null;
    }

    public double getData(Identifier taxon, Trait trait) {
        return getData(whichTaxon(taxon), whichTrait(trait));
    }

    public double getData(int taxon, int trait) {
        int pheno = phenotypeColumn[trait][0];
        if (phenotypeRow[taxon][pheno] == -1) {
            return Double.NaN;
        } else {
            return phenotypes[pheno].getData(phenotypeRow[taxon][pheno], phenotypeColumn[trait][1]);
        }
    }

    public void setData(Identifier taxon, Trait trait, double value) {
        setData(whichTaxon(taxon), whichTrait(trait), value);
    }

    public void setData(int taxon, int trait, double value) {
        int pheno = phenotypeColumn[trait][0];
        phenotypes[pheno].setData(phenotypeRow[taxon][pheno], phenotypeColumn[trait][1], value);
    }

    public SimplePhenotype simpleCopy() {
        return new SimplePhenotype(getTaxa(), getTraits(), getData());
    }
}
