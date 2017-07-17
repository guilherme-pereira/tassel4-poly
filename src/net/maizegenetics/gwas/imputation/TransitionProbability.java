package net.maizegenetics.gwas.imputation;

public class TransitionProbability {
	protected double[][] probabilityOfATransition;
	protected double[][] adjustedProbability;
	protected int[] positions;
	double avgSegmentLength = Double.NaN;
	
	/**
	 * Use this function to set the transition probability matrix when it is to be used without specifying nodes for getTransitionProbability
	 * @param probabilityMatrix the matrix of transition probabilities, rows = state1, columns = state2
	 */
	public void setTransitionProbability(double[][] probabilityMatrix) {
		probabilityOfATransition = probabilityMatrix;
	}
	
	public double getTransitionProbability(int state1, int state2) {
		return adjustedProbability[state1][state2];
//		return probabilityOfATransition[state1][state2];
	}
	
	public double getLnTransitionProbability(int state1, int state2) {
		return Math.log(getTransitionProbability(state1, state2));
	}
	
	public int getNumberOfStates() {
		return probabilityOfATransition.length;
	}
	
	/**
	 * @param transitionCounts the transition counts for this set of 
	 * @param chromosomeLength
	 * @param numberOfTaxa
	 */
	public void setTransitionCounts(int[][] transitionCounts, int chromosomeLength, int numberOfTaxa) {
		int n = transitionCounts.length;
		int totalTransitions = 0;
		int[] rowSums = new int[n];
		int rowCount = 0;
		for (int[] row:transitionCounts) {
			for (int cell: row) {
				totalTransitions += cell;
				rowSums[rowCount] += cell;
			}
			rowCount++;
		}
		
		avgSegmentLength = ((double) chromosomeLength) * ((double) numberOfTaxa) / ((double) totalTransitions);
		
		probabilityOfATransition = new double[n][n];
		for (int row = 0; row < n; row++) {
			for (int col = 0; col < n; col++) {
				probabilityOfATransition[row][col] = ((double) transitionCounts[row][col]) / ((double) rowSums[row]);
			}
		}
	}

	public void setNode(int node) {
		if (node <= 0) return;
		int n = probabilityOfATransition.length;
		adjustedProbability = new double[n][n];
		int segmentLength = Math.abs(positions[node] - positions[node - 1]);
		double m;
		for (int row = 0; row < n; row++) {
			double offdiagsum = 0;
			for (int col = 0; col < n; col++) {
				if (col != row) {
					m = -Math.log(1 - 2 * probabilityOfATransition[row][col]) * segmentLength / avgSegmentLength / 2;
					adjustedProbability[row][col] = (1 - Math.exp(-2*m)) / 2;
					offdiagsum += adjustedProbability[row][col];
				}
			}
			adjustedProbability[row][row] = 1 - offdiagsum;
		}
	}

	public void setAverageSegmentLength(double length) { avgSegmentLength = length; }
	public void setPositions(int[] positions) { this.positions = positions; }
}
