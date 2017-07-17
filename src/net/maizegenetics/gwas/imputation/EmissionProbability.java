package net.maizegenetics.gwas.imputation;

public class EmissionProbability {
	double[][] probObsGivenState; //state in rows, observation in columns
	
	public EmissionProbability() {
		
	}
	
	public double getProbObsGivenState(int state, int obs) {
		return probObsGivenState[state][obs];
	}
	
	public double getLnProbObsGivenState(int state, int obs) {
		return Math.log(getProbObsGivenState(state, obs));
	}

	public double getProbObsGivenState(int state, int obs, int node) {
		return probObsGivenState[state][obs];
	}
	
	public double getLnProbObsGivenState(int state, int obs, int node) {
		return Math.log(getProbObsGivenState(state, obs, node));
	}

	public void setEmissionProbability(double[][] probabilityMatrix) {
		probObsGivenState = probabilityMatrix;
	}
}
