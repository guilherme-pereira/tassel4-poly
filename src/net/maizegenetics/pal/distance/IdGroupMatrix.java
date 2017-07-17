package net.maizegenetics.pal.distance;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.report.TableReport;

import java.io.Serializable;

/**
 * Created by IntelliJ IDEA.
 * User: ed
 * Date: Jan 10, 2007
 * Time: 12:13:28 PM
 * To change this template use File | Settings | File Templates.
 */
public interface IdGroupMatrix extends Serializable, IdGroup, TableReport {

    /** returns representation of this alignment as a string */
    String toString();

    	/**
	 * Returns the number of rows and columns that the distance matrix has.
	 */
    int getSize();

    	/**
	 * Returns the distances as a 2-dimensional array of doubles. Matrix is cloned first so it can be altered freely.
	 */
    double[][] getClonedDistances();

    double getDistance(int row, int col);

    	/**
	 * Returns the mean pairwise distance of this matrix
	 */
    double meanDistance();

    	/**
	 * test whether this matrix is a symmetric distance matrix
	 *
	 */
    boolean isSymmetric();
}
