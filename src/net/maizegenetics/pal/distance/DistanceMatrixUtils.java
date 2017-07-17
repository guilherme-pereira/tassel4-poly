// DistanceMatrixUtils.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.distance;

import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;

import java.io.Serializable;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;


/**
 * Auxillary functions for distance matrices<p>
 *
 * @version $Id: DistanceMatrixUtils.java,v 1.1 2007/01/12 03:26:14 tcasstevens Exp $
 *
 * @author Alexei Drummond
 */
public class DistanceMatrixUtils implements Serializable {

	/**
	 * compute squared distance to second distance matrix.
	 * If both matrices have the same size it is assumed that the order of the taxa
	 * is identical.
	 */
	public static double squaredDistance(DistanceMatrix mat1, DistanceMatrix mat2, boolean weighted) {

		boolean aliasNeeded = false;
		if (mat1.getSize() != mat2.getSize())
		{
			aliasNeeded = true;
		}

		int[] alias = null;

		if (aliasNeeded) {
			if (mat1.getSize() > mat2.getSize()) {
				//swap so mat1 is the smaller of the two
				DistanceMatrix temp = mat2;
				mat2 = mat1;
				mat1 = temp;
			}
			alias = new int[mat1.getSize()];
			for (int i = 0; i < alias.length; i++) {
				alias[i] = mat2.whichIdNumber(mat1.getIdentifier(i).getName());
			}
		} else {
			alias = new int[mat1.getSize()];
			for (int i = 0; i < alias.length; i++) {
				alias[i] = i;
			}
		}


		double sum = 0;
		int ai;
		final double[][] mat1Distance = mat1.getDistances();
		final double[][] mat2Distance = mat2.getDistances();
		for (int i = 0; i < mat1.getSize()-1; i++)
		{
			ai = alias[i];

			for (int j = i+1; j < mat1.getSize(); j++)
			{
				double diff = mat1Distance[i][j] - mat2Distance[ai][alias[j]];
				double weight;
				if (weighted)
				{
					// Fitch-Margoliash weight
					// (variances proportional to distances)
					weight = 1.0/(mat1Distance[i][j]*mat2Distance[ai][alias[j]]);
				}
				else
				{
					// Cavalli-Sforza-Edwards weight
					// (homogeneity of variances)
					weight = 1.0;
				}
				sum += weight*diff*diff;
			}
		}

		return 2.0*sum; // we counted only half the matrix
	}

	/**
	 * Returns a distance matrix with the specified taxa removed.
	 */
	public static DistanceMatrix minus(DistanceMatrix parent, int taxaToRemove) {

		int size = parent.getIdCount() - 1;

		double[][] distances = new double[size][size];
		Identifier[] ids = new Identifier[size];
		int counti = 0, countj = 0;
		for (int i = 0; i < size; i++) {
			if (counti == taxaToRemove) {
				counti += 1;
			}
			ids[i] = parent.getIdentifier(counti);

			countj = 0;
			final double[][] parentDistance = parent.getDistances();
			for (int j = 0; j < size; j++) {
				if (countj == taxaToRemove) {
					countj += 1;
				}
				distances[i][j] = parentDistance[counti][countj];
				countj += 1;
			}
			counti += 1;
		}

		DistanceMatrix smaller = new DistanceMatrix(distances, new SimpleIdGroup(ids));

		return smaller;
	}
        
    /**
     * Calculates the IBS distance between two taxa with bitsets for for major and minor allele
     * @param iMajor
     * @param iMinor
     * @param jMajor
     * @param jMinor
     * @return 
     */ 
    public static double getIBSDistance(long[] iMajor, long[] iMinor, long[] jMajor, long[] jMinor) {
        int sameCnt=0, diffCnt=0, hetCnt=0;
        for(int x=0; x<iMajor.length; x++) {
            long same=(iMajor[x]&jMajor[x])|(iMinor[x]&jMinor[x]);
            long diff=(iMajor[x]&jMinor[x])|(iMinor[x]&jMajor[x]);
            long hets=same&diff;
            sameCnt+=BitUtil.pop(same);
            diffCnt+=BitUtil.pop(diff);
            hetCnt+=BitUtil.pop(hets);
        }
        double identity=(double)(sameCnt+(hetCnt/2))/(double)(sameCnt+diffCnt+hetCnt);
        double dist=1-identity;
        return dist;
    }
    
    public static double getIBSDistance(BitSet iMajor, BitSet iMinor, BitSet jMajor, BitSet jMinor) {
        return getIBSDistance(iMajor.getBits(), iMinor.getBits(), jMajor.getBits(), jMinor.getBits());
    }

}
