package net.maizegenetics.gwas.imputation;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;

public class PhasedEmissionProbability extends EmissionProbability {
	float[][][] parentHaplotypeProbability; //first dimension is the site, second dimension is the haplotype (0-3), third dimension is nucleotide (A=0, C=1, G=2, T=3)
	
	public PhasedEmissionProbability() {

		
	}
	
	@Override
	public double getProbObsGivenState(int state, int obs, int node) {
		//state is an ordered pair of haplotypes, obs is the observed genotype, node is the site
		//P(obs=X|state=H1H2) = P(obs=X|state=H1H2, H1=A, H2=A)P(H1=A)P(H2=A) + P(obs=X|state=H1H2, H1=A, H2=C)P(H1=A)P(H2=C) + ... + P(obs=X|state=H1H2, H1=T, H2=T)P(H1=T)P(H2=T)
		double prob = 0.0;
		byte[] haplotypes = AlignmentUtils.getDiploidValues((byte) state);
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				prob += conditionalProbability(obs, i, j) * parentHaplotypeProbability[node][haplotypes[0]][i] * parentHaplotypeProbability[node][haplotypes[1]][j];
			}
		}
		return prob;
	}

	@Override
	public double getLnProbObsGivenState(int state, int obs, int node) {
		return Math.log(getProbObsGivenState(state, obs, node));
	}

	private double conditionalProbability(int obs, int N1, int N2) {
		double p = 0.0;
		
		return p;
	}
	
	public void setParentHaplotypeProbability(float[][][] probability) {
		parentHaplotypeProbability = probability;
	}
	
	public void setParentHaplotypeProbability(Alignment a) {
		
	}
	
	public void setParentHaplotypeProbability(PopulationData pop) {
		
	}

	
}
