// DistanceMatrix.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.distance;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.io.FormattedOutput;
import net.maizegenetics.pal.report.TableReport;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;


/**
 * storage for pairwise distance matrices.<p>
 *
 * features:
 * - printing in in PHYLIP format,
 * - computation of (weighted) squared distance to other distance matrix
 * - Fills in all of array...
 *
 * @version $Id: DistanceMatrix.java,v 1.8 2009/07/07 16:19:37 tcasstevens Exp $
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */
public class DistanceMatrix implements IdGroupMatrix, TableReport {

    //
    // Private stuff
    //
    /** sequence identifiers */
    private IdGroup idGroup;
    /** distances [seq1][seq2] */
    private double[][] distance = null;
    static final long serialVersionUID = 4725925229860707633L;

    /** I like doing things my self! */
    private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
        out.writeByte(1); //Version number
        out.writeObject(idGroup);
        out.writeObject(distance);
    }

    private void readObject(java.io.ObjectInputStream in) throws IOException, ClassNotFoundException {
        byte version = in.readByte();
        switch (version) {
            default: {
                idGroup = (IdGroup) in.readObject();
                distance = (double[][]) in.readObject();
            }
        }
    }

    /** constructor */
    public DistanceMatrix() {
    }

    /** constructor taking distances array and IdGroup */
    public DistanceMatrix(double[][] distance, IdGroup idGroup) {
        super();
        this.distance = distance;
        this.idGroup = idGroup;
    }

    /**
     * constructor that takes a distance matrix and clones the distances
     * and IdGroup
     */
    public DistanceMatrix(DistanceMatrix dm) {
        distance = net.maizegenetics.pal.util.Utils.getCopy(dm.distance);
        idGroup = SimpleIdGroup.getInstance(dm.idGroup);
    }

    /**
     * constructor that takes a distance matrix and clones the distances,
     * of a the identifiers in idGroup.
     */
    public DistanceMatrix(DistanceMatrix dm, IdGroup subset) {

        int index1, index2;

        distance = new double[subset.getIdCount()][subset.getIdCount()];
        for (int i = 0; i < distance.length; i++) {
            index1 = dm.whichIdNumber(subset.getIdentifier(i).getName());
            distance[i][i] = dm.distance[index1][index1];
            for (int j = 0; j < i; j++) {
                index2 = dm.whichIdNumber(subset.getIdentifier(j).getName());
                distance[i][j] = dm.distance[index1][index2];
                distance[j][i] = distance[i][j];
            }
        }
        idGroup = subset;
    }

    /** print alignment (PHYLIP format) */
    public void printPHYLIP(PrintWriter out) {
        // PHYLIP header line
        out.println("  " + distance.length);
        FormattedOutput format = FormattedOutput.getInstance();

        for (int i = 0; i < distance.length; i++) {
            format.displayLabel(out,
                    idGroup.getIdentifier(i).getName(), 10);
            out.print("      ");

            for (int j = 0; j < distance.length; j++) {
                // Chunks of 6 blocks each
                if (j % 6 == 0 && j != 0) {
                    out.println();
                    out.print("                ");
                }

                out.print("  ");
                format.displayDecimal(out, distance[i][j], 5);
            }
            out.println();
        }
    }

    /** returns representation of this alignment as a string */
    public String toString() {

        StringWriter sw = new StringWriter();
        printPHYLIP(new PrintWriter(sw));

        return sw.toString();
    }

    /** compute squared distance to second distance matrix */
    public double squaredDistance(DistanceMatrix mat, boolean weighted) {
        double sum = 0;
        for (int i = 0; i < distance.length - 1; i++) {
            for (int j = i + 1; j < distance.length; j++) {
                double diff = distance[i][j] - mat.distance[i][j];
                double weight;
                if (weighted) {
                    // Fitch-Margoliash weight
                    // (variances proportional to distances)
                    weight = 1.0 / (distance[i][j] * distance[i][j]);
                } else {
                    // Cavalli-Sforza-Edwards weight
                    // (homogeneity of variances)
                    weight = 1.0;
                }
                sum += weight * diff * diff;
            }
        }

        return 2.0 * sum; // we counted only half the matrix
    }

    /** compute absolute distance to second distance matrix */
    public double absoluteDistance(DistanceMatrix mat) {
        double sum = 0;
        for (int i = 0; i < distance.length - 1; i++) {
            for (int j = i + 1; j < distance.length; j++) {
                double diff =
                        Math.abs(distance[i][j] - mat.distance[i][j]);

                sum += diff;
            }
        }

        return 2.0 * sum; // we counted only half the matrix
    }

    /**
     * Returns the number of rows and columns that the distance matrix has.
     */
    public int getSize() {
        return distance.length;
    }

    /**
     * Returns the distances as a 2-dimensional array of doubles. Matrix is cloned first so it can be altered freely.
     */
    public final double[][] getClonedDistances() {
        return net.maizegenetics.pal.util.Utils.getCopy(distance);
    }

    /**
     * Returns the distances as a 2-dimensional array of doubles (in the actual array used to store the distances)
     */
    protected final double[][] getDistances() {
        return net.maizegenetics.pal.util.Utils.getCopy(distance);
    }

    public final double getDistance(final int row, final int col) {
        return distance[row][col];
    }

    /**
     * Sets both upper and lower triangles.
     */
    public void setDistance(int i, int j, double dist) {
        distance[i][j] = distance[j][i] = dist;
    }

    /**
     * Adds a delta to both upper and lower triangle distances.
     */
    public void addDistance(int i, int j, double delta) {
        distance[i][j] += delta;
        distance[j][i] += delta;
    }

    /**
     * Returns the mean pairwise distance of this matrix
     */
    public double meanDistance() {
        double dist = 0.0;
        int count = 0;
        for (int i = 0; i < distance.length; i++) {
            for (int j = 0; j < distance[i].length; j++) {
                if ((i != j)&&(!Double.isNaN(distance[i][j]))) {
                    dist += distance[i][j];
                    count += 1;
                }
            }
        }
        return dist / (double) count;
    }

    //IdGroup interface
    public Identifier getIdentifier(int i) {
        return idGroup.getIdentifier(i);
    }

    public void setIdentifier(int i, Identifier ident) {
        idGroup.setIdentifier(i, ident);
    }

    public int getIdCount() {
        return idGroup.getIdCount();
    }

    public int whichIdNumber(String name) {
        return idGroup.whichIdNumber(name);
    }

    public int whichIdNumber(Identifier id) {
        return idGroup.whichIdNumber(id);
    }

    /**
     * Return id group of this alignment.
     * @deprecated distance matrix now implements IdGroup
     */
    public IdGroup getIdGroup() {
        return idGroup;
    }

    /**
     * test whether this matrix is a symmetric distance matrix
     *
     */
    public boolean isSymmetric() {
        for (int i = 0; i < distance.length; i++) {
            if (distance[i][i] != 0) {
                return false;
            }
        }
        for (int i = 0; i < distance.length - 1; i++) {
            for (int j = i + 1; j < distance.length; j++) {
                if (distance[i][j] != distance[j][i]) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @param fromID the thing (taxa,sequence) from which we want to find the closest (excluding self)
     * @param exlcusion indexes of things that should not be considered, may be null
     * @return the index of the thing closest to the specified
     * @note if fromID not a valid name then return -1
     */
    public int getClosestIndex(String fromID, String[] exclusion) {
        int index = whichIdNumber(fromID);
        if (index < 0) {
            return -1;
        }
        int[] exclusionIndexes;
        if (exclusion == null) {
            exclusionIndexes = null;
        } else {
            exclusionIndexes = new int[exclusion.length];
            for (int i = 0; i < exclusion.length; i++) {
                exclusionIndexes[i] = whichIdNumber(exclusion[i]);
            }
        }
        return getClosestIndex(index, exclusionIndexes);
    }

    private final boolean isIn(int value, int[] set) {
        if (set == null) {
            return false;
        }
        for (int i = 0; i < set.length; i++) {
            if (set[i] == value) {
                return true;
            }
        }
        return false;
    }

    /**
     * @param fromIndex the index of the thing (taxa,sequence) from which we want to find the closest (excluding self)
     * @param exlcusion indexes of things that should not be considered, may be null
     * @return the index of the member closes to the specified
     */
    public int getClosestIndex(int fromIndex, int[] exclusion) {
        double min = Double.POSITIVE_INFINITY;
        int index = -1;
        for (int i = 0; i < distance.length; i++) {
            if (i != fromIndex && !isIn(i, exclusion)) {
                double d = distance[fromIndex][i];
                if (d < min) {
                    min = d;
                    index = i;
                }
            }
        }
        return index;
    }

    protected final void setIdGroup(IdGroup base) {
        this.idGroup = SimpleIdGroup.getInstance(base);
    }

    protected final void setDistances(double[][] matrix) {
        this.distance = matrix;
    }

    public Object[] getTableColumnNames() {
        String[] colNames = new String[getSize() + 1];
        colNames[0] = "Taxa";
        for (int i = 0; i < distance[0].length; i++) {
            colNames[i + 1] = getIdentifier(i).toString();
        }
        return colNames;
    }

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    public Object[] getRow(int row) {

        Object[] result = new Object[distance[row].length + 1];
        result[0] = getIdentifier(row);
        for (int j = 1; j <= distance[row].length; j++) {
            result[j] = "" + distance[row][j - 1];
        }

        return result;

    }

    public String getTableTitle() {
        return "Alignment Distance Matrix";
    }

    public int getRowCount() {
        if (distance != null) {
            return distance.length;
        } else {
            return 0;
        }
    }

    public int getElementCount() {
        return getRowCount() * getColumnCount();
    }

    public int getColumnCount() {
        if ((distance != null) && distance[0] != null) {
            return distance[0].length + 1;
        } else {
            return 0;
        }
    }

	@Override
	public Object[][] getTableData() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Object[][] getTableData(int start, int end) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		if (columnIndex == 0) return getIdentifier(rowIndex);
		return new Double(distance[rowIndex][columnIndex - 1]);
	}

    public String getColumnName(int col) {
        if (col == 0) {
            return "Taxa";
        }
        return getIdentifier(col-1).toString();
    }

}
