package net.maizegenetics.gwas.NAM;

import java.util.ArrayList;

import net.maizegenetics.jGLiM.LinearModelUtils;
import net.maizegenetics.jGLiM.dm.CovariateModelEffect;
import net.maizegenetics.jGLiM.dm.LinearModelForStepwiseRegression;
import net.maizegenetics.jGLiM.dm.ModelEffect;
import net.maizegenetics.jGLiM.dm.SweepFastLinearModel;

public class ForwardRegressionModel {
	ArrayList<ModelEffect> modelEffects = new ArrayList<ModelEffect>();
	LinearModelForStepwiseRegression lmsr;
	double bestModelSS = 0;
	SnpInfo bestsnp = null;
	double enterlimit;
	int maxsnps;
	int nsnpsInModel = 0;
	boolean wasUpdated = true;
	
	public ForwardRegressionModel(ArrayList<ModelEffect> initialEffects, double[] phenotype, double enterlimit, int maxsnps) {
		modelEffects.addAll(initialEffects);
		this.enterlimit = enterlimit;
		this.maxsnps = maxsnps;
		lmsr = new LinearModelForStepwiseRegression(modelEffects, phenotype);
	}
	
	/**
	 * Call this function to test another snp.
	 * @param snp	the snp info for this snp
	 */
	public void addNextSnp(SnpInfo snp) {
		double modelss = lmsr.testNewEffect(snp.genotype);
		if (modelss > bestModelSS) {
			bestModelSS = modelss;
			bestsnp = snp;
		}
	}
	
	/**
	 * Call this function when all the snps have been tested once. 
	 * If the p-value for the most significant SNP is less than the entry p-value and
	 * the maximum number of snps has not been exceeded,
	 * then the SNP will be added to the model and the function will return true 
	 * @return	true if a SNP was added to the model, false otherwise.
	 */
	public boolean updateModel() {
		double[] Fp = lmsr.getFpFromModelSS(bestModelSS);
		if (Fp[1] < enterlimit) {
			lmsr.addEffect(new CovariateModelEffect(bestsnp.genotype, bestsnp));
			nsnpsInModel++;
			System.out.println(this.toString());
			System.out.println();
			bestsnp = null;
			bestModelSS = 0;
			wasUpdated = true;
			if (nsnpsInModel >= maxsnps) {
				wasUpdated = false;
				return false;
			}
			wasUpdated = true;
			return true;
		} else {
			wasUpdated = false;
			return false;
		}
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		String sep = "\t";
		SweepFastLinearModel sflm = lmsr.getLinearModel();
		double[] beta = sflm.getBeta();
		int neffects = lmsr.getModelEffects().size();
		int nBeta = beta.length;
		int start = nBeta - neffects;
		double[] errorSSdf = sflm.getResidualSSdf();
		for (int i = 2; i < neffects; i++) {
			if (i > 2) sb.append("\n");
			SnpInfo id = (SnpInfo) modelEffects.get(i).getID();
			sb.append(id.chromosome);
			sb.append(sep).append(id.pos);
			sb.append(sep).append(id.allele);
			sb.append(sep).append(beta[i + start]);

			//calculate f and p
			double[] ssdf = sflm.getMarginalSSdf(i);
			double F = ssdf[0] / ssdf[1] / errorSSdf[0] * errorSSdf[1];
			double p;
			try {
				p = LinearModelUtils.Ftest(F, ssdf[1], errorSSdf[1]);
			} catch (Exception e) {
				p = Double.NaN;
			}

			sb.append(sep).append(F);
			sb.append(sep).append(p);
			sb.append(sep).append(-Math.log10(p));
		}
		return sb.toString();
	}
	
	public boolean wasUpdated() {return wasUpdated;}
}
