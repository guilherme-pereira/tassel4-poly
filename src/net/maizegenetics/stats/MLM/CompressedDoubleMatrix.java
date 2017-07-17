package net.maizegenetics.stats.MLM;

import java.util.ArrayList;
import java.util.Collections;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.pal.tree.Tree;
import net.maizegenetics.pal.tree.TreeClusters;
import net.maizegenetics.pal.tree.UPGMATree;

public class CompressedDoubleMatrix {
	protected DistanceMatrix kinshipMatrix;
	protected TreeClusters theClusters;
	protected int numberOfGroups;
	protected DoubleMatrix adjustmentMatrix;
	protected DoubleMatrix compressedMatrix;
	protected int[] distance2clusterIndex;
	protected int[][] groupMembers;
	Tree myTree;
	
	public enum kinshipMethod {avg, min, max, median};
	
	public CompressedDoubleMatrix(DistanceMatrix kinship, Tree theTree) {
		kinshipMatrix = kinship;
		myTree = theTree;
		init();
	}
	
	protected void init() {
		theClusters = new TreeClusters(myTree);
		
		int nrow = kinshipMatrix.getSize();
		
		distance2clusterIndex = new int[nrow];
		for (int i = 0; i < nrow; i++) distance2clusterIndex[i] = myTree.whichIdNumber(kinshipMatrix.getIdentifier(i));
	}

	public int getNumberOfGroups() {
		return numberOfGroups;
	}

	/**
	 * This sets the number of groups and calculates the output matrices. 
	 * It must be run to set the number of groups before retrieving the compressed and adjustment matrices. 
	 * @param numberOfGroups the number of groups into which the taxa should be compressed
	 * 	 
	 * */
	public void setNumberOfGroups(int numberOfGroups) {
		int nrow = kinshipMatrix.getSize();

		int[] clusterIndex = theClusters.getGroups(numberOfGroups);
		
		int maxIndex = 0;
		for (int i : clusterIndex) maxIndex = Math.max(maxIndex, i);
		int ncol = maxIndex + 1;
		this.numberOfGroups = ncol;
		
		adjustmentMatrix = DoubleMatrixFactory.DEFAULT.make(nrow, ncol, 0);
		for (int i = 0; i < nrow; i++) {
			int thiscol = clusterIndex[distance2clusterIndex[i]];
			adjustmentMatrix.set(i, thiscol, 1);
		}
		
		groupMembers = new int[ncol][];
		for (int i = 0; i < ncol; i++) {
			int n = (int) adjustmentMatrix.columnSum(i);
			groupMembers[i] = new int[n];
			int count = 0;
			for (int r = 0; r < nrow; r++) {
				if (adjustmentMatrix.get(r, i) == 1) groupMembers[i][count++] = r; 
			}
		}
		
	}
	
	/**
	 * @param km the kinship method used to calculate the group kinships
	 * @return the compressed distance matrix for the set number of groups
	 */
	public DoubleMatrix getCompressedMatrix(kinshipMethod km) {
		int ngroups = adjustmentMatrix.numberOfColumns();
		compressedMatrix = DoubleMatrixFactory.DEFAULT.make(ngroups, ngroups, 0);


		for (int g1 = 0; g1 < ngroups; g1++) {
			for (int g2 = g1; g2 < ngroups; g2++) {
				ArrayList<Double> kinshipCoefficients = new ArrayList<Double>();
				for (int t1: groupMembers[g1]) {
					for (int t2: groupMembers[g2]) {
						kinshipCoefficients.add(kinshipMatrix.getDistance(t1, t2));
					}
				}
				
				double groupDistance = 0;
				switch(km) {
				case avg:
					double n = kinshipCoefficients.size();
					double total = 0;
					for (Double kc:kinshipCoefficients) total += kc;
					groupDistance = total / n;
					break;
				case min:
					groupDistance = kinshipCoefficients.get(0);
					for (Double kc:kinshipCoefficients) groupDistance = Math.min(groupDistance, kc);
					break;
				case max:
					groupDistance = kinshipCoefficients.get(0);
					for (Double kc:kinshipCoefficients) groupDistance = Math.max(groupDistance, kc);
					break;
				case median:
					Collections.sort(kinshipCoefficients);
					int ncoeff = kinshipCoefficients.size();
					
					if (ncoeff % 2 == 0) { //ncoeff is even
						int midpoint = ncoeff / 2;
						groupDistance = (kinshipCoefficients.get(midpoint) + kinshipCoefficients.get(midpoint - 1))/2;
					} else { //ncoeff is odd
						int midpoint = (ncoeff - 1)/2;
						groupDistance = kinshipCoefficients.get(midpoint);
					}
					break;
				}
				compressedMatrix.set(g1, g2, groupDistance);
				compressedMatrix.set(g2, g1, groupDistance);
			}
		}
		return compressedMatrix;
	}
	
	/**
	 * @return the matrix used to post-multiply the Z matrix to produce the compressed Z matrix
	 */
	public DoubleMatrix getAdjustmentMatrix() {
		return adjustmentMatrix;
	}
	
	public DoubleMatrix getCompressedZKZ(DoubleMatrix Z, kinshipMethod km) {
		DoubleMatrix compressedZ = getCompressedZ(Z);
		getCompressedMatrix(km);
		return compressedZ.mult(compressedMatrix).tcrossproduct(compressedZ);
	}
	
	public DoubleMatrix getCompressedZ(DoubleMatrix originalZ) {
		return originalZ.mult(adjustmentMatrix);
	}
	
}
