package net.maizegenetics.gwas.imputation;

public class ViterbiAlgorithm {

	//adapted from Rabiner Proceedings of the IEEE 77(2):257-286
	//initialize
	//d(0, i) = p(Obs-0|S(i)*p(S(i))
	//where d0 = path length, Obs-0 = observation 0, S0 = true state 0
	//iterate:
	//for t = 1 to n
	//d(t,S(j)) = max(j){d(t-1,S(i)) * p[S(j)|S(i)] * p[Obs(t)|S(j)]}
	//h(t, j) = the value of i that maximizes distance
	//where h() = path history
	//termination
	//choose state that maximizes path length
	//back tracking
	//S(t) = h(t+1, S(t+1)), to decode best sequence
	
	TransitionProbability myTransitionMatrix;
	EmissionProbability probObservationGivenState;
	byte[] obs;
	byte[][] history;
	double[] distance;
	int numberOfStates; //number of true states;
	double[] probTrueStates; //ln of probabilities
	int numberOfObs;
	byte[] finalState;
	
	
	public ViterbiAlgorithm(byte[] observations, TransitionProbability transitionMatrix, EmissionProbability obsGivenTrue, double[] pTrue) {
		obs = observations;
		numberOfObs = obs.length;
		numberOfStates = transitionMatrix.getNumberOfStates();
		
		myTransitionMatrix = transitionMatrix;
		probObservationGivenState = obsGivenTrue;
		probTrueStates = new double[pTrue.length];
		for (int i = 0; i < pTrue.length; i++) {
			probTrueStates[i] = Math.log(pTrue[i]);
		}
		
		history = new byte[numberOfStates][numberOfObs];
		distance = new double[numberOfStates];
	}
	
	public void calculate() {
		initialize();
		for (int i = 1; i < numberOfObs; i++) {
			updateDistanceAndHistory(i);
		}
	}
	
	public void initialize() {
		for (int i = 0; i < numberOfStates; i++) {
			try{
				distance[i] = probObservationGivenState.getLnProbObsGivenState(i, obs[0], 0) + probTrueStates[i];
			} catch(Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	public void updateDistanceAndHistory(int node) {
		double[][] candidateDistance = new double[numberOfStates][numberOfStates];
		myTransitionMatrix.setNode(node);
		for (int i = 0; i < numberOfStates; i++) {
			for (int j = 0; j < numberOfStates; j++) {
				candidateDistance[i][j] = distance[i] + myTransitionMatrix.getLnTransitionProbability(i, j) + probObservationGivenState.getLnProbObsGivenState(j, obs[node], node);
			}
		}

		//find the maxima
		int[] max = new int[numberOfStates];
		for (int i = 0; i < numberOfStates; i++) {
			for (int j = 0; j < numberOfStates; j++) {
				if (candidateDistance[i][j] > candidateDistance[max[j]][j]) max[j] = i;
			}
		}

		//update distance and history
		for (int j = 0; j < numberOfStates; j++) {
			distance[j] = candidateDistance[max[j]][j];
			history[j][node] = (byte) max[j];
		}
		
		//if the min distance is less than -1e100, subtract the max distance;
		double maxd = distance[0];
		double mind = 0;
		for (int i = 0; i < numberOfStates; i++) {
			if (distance[i] > maxd) maxd = distance[i];
			if (distance[i] != Double.NEGATIVE_INFINITY && distance[i] < mind) mind = distance[i];
		}
		if (mind < -1e100) {
			for (int i = 0; i < numberOfStates; i++) {
				distance[i] -= maxd;
			}
		}
	}
	
	//decode the most probable state sequence
	public byte[] getMostProbableStateSequence() {
		byte[] seq = new byte[numberOfObs];
		byte finalState = 0;
		for (int i = 1; i < numberOfStates; i++) {
			if (distance[i] > distance[finalState]) finalState = (byte) i;
		}
		
		//S(t) = h(t+1, S(t+1)), to decode best sequence
		seq[numberOfObs - 1] = finalState;
		for (int i = numberOfObs - 2; i >= 0; i--) {
			 seq[i] = history[seq[i + 1]][i + 1];
		}
		return seq;
	}
		
	public void setStateProbability(double[] probTrueState) {
		int n = probTrueState.length;
		probTrueStates = new double[n];
		for (int i = 0; i < n; i++) {
			probTrueStates[i] = Math.log(probTrueState[i]);
		}
	}
}
