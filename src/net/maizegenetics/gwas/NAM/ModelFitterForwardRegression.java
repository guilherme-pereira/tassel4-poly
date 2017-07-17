 package net.maizegenetics.gwas.NAM;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import net.maizegenetics.jGLiM.LinearModelUtils;
import net.maizegenetics.jGLiM.dm.FactorModelEffect;
import net.maizegenetics.jGLiM.dm.LinearModelForStepwiseRegression;
import net.maizegenetics.jGLiM.dm.ModelEffect;
import net.maizegenetics.jGLiM.dm.CovariateModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffectUtils;
import net.maizegenetics.jGLiM.dm.SweepFastLinearModel;
	
public class ModelFitterForwardRegression extends ModelFitter {
	double[] data;
	double enterLimit = 1e-4;
	LinearModelForStepwiseRegression lmsr;
	
	public ModelFitterForwardRegression(int chromosome, FileNames files) {
		super(chromosome, files);
		if (!Double.isNaN(files.enterlimit)) enterLimit = files.enterlimit;
	}

	@Override
	protected void testSnps() {
		initializeStepFile();

		//get data for this chromosome
		int chromosome = files.chromosome;
		int nSamples = residuals.length;
		data = new double[nSamples];
		for (int i = 0; i < nSamples; i++) data[i] = residuals[i][chromosome - 1];
		
		//create the base model
		lmsr = getBaseModel();
		SnpInfo nextSnp = findNextTerm();
		while (nextSnp.p < enterLimit) {
			writeSNPToStepFile(nextSnp);
			lmsr.addEffect(new CovariateModelEffect(nextSnp.genotype, nextSnp));
			nextSnp = findNextTerm();
		}
		
		writeOutputToFile();
	}
	
	//for the next term (max modelss) need pos, allele, projection, F, p
	protected SnpInfo findNextTerm() {
		int chromosome = files.chromosome;
		SnpInfo bestSnp = null;
		double bestSS = 0;
		snpdata.reset();
		while (snpdata.next()) {
			double[] parents = snpdata.getGenotype();
			int pos = snpdata.getPosition();
			
			//build and solve the model
			double[] snp = projectSnp(parents, pos);
			double ms =  lmsr.testNewEffect(snp);
			
			if (ms > bestSS) {
				bestSS = ms;
				double[] Fp = lmsr.getFpFromModelSS(ms);
				bestSnp = new SnpInfo(chromosome, pos, snpdata.getAllele(), snp, Fp[0], Fp[1]);
			}
		}
		System.out.println(chromosome + ", " + bestSnp.pos + ", " + bestSnp.allele + ", " + bestSnp.F + ", " + bestSnp.p);
		return bestSnp;
	}
	
	protected LinearModelForStepwiseRegression getBaseModel() {
		int nSamples = residualSamples.length;
		int[] mean = new int[nSamples];
		
		//set up the mean and pop model effects
		ArrayList<ModelEffect> effects = new ArrayList<ModelEffect>();
		ModelEffect memean = new FactorModelEffect(mean, false);
		effects.add(memean);
		ModelEffect mepop = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(pops), true);
		effects.add(mepop);
		
		//create the base model
		return new LinearModelForStepwiseRegression(effects, data);
	}
	
	public void writeOutputToFile() {
		//calculate fit model information
		SweepFastLinearModel sflm = lmsr.getLinearModel();
		double[] beta = sflm.getBeta();
		ArrayList<ModelEffect> modelEffects = lmsr.getModelEffects();
		int neffects = modelEffects.size();
		double[] errorSSdf = sflm.getResidualSSdf();
		BufferedWriter bw = null;
		
		try {
			bw = new BufferedWriter(new FileWriter(files.chrmodel));
			bw.write("Chr\tposition\tallele\teffect\tF\tp\tlog10p");
			bw.newLine();
		} catch (IOException e1) {
			e1.printStackTrace();
			System.exit(-1);
		}
		String sep = "\t";
		
		
		int nBeta = beta.length;
		int start = nBeta - neffects;
		int chromosome = files.chromosome;
		for (int i = 2; i < neffects; i++) {
			StringBuilder sb = new StringBuilder();
			SnpInfo id = (SnpInfo) modelEffects.get(i).getID();
			sb.append(chromosome);
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
			
			try {
				bw.write(sb.toString());
				bw.newLine();
			} catch (IOException e1) {
				System.out.println("File write failed at:\n");
//				System.out.println(sb.toString());
			}
			System.out.println(sb.toString());
			
		}
		
		try {
			bw.close();
		} catch (IOException e1) {
			e1.printStackTrace();
			System.exit(-1);
		}
		
	}
	
	public void initializeStepFile() {
		BufferedWriter bw;
		try {
			bw = new BufferedWriter(new FileWriter(files.chrsteps));
			bw.write("Chromosome\tposition\tallele\tF\tp");
			bw.newLine();
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public void writeSNPToStepFile(SnpInfo snp) {
		BufferedWriter bw;
		try {
			bw = new BufferedWriter(new FileWriter(files.chrsteps, true));
			bw.write(Integer.toString(snp.chromosome));
			bw.write("\t");
			bw.write(Integer.toString(snp.pos));
			bw.write("\t");
			bw.write(snp.allele);
			bw.write("\t");
			bw.write(Double.toString(snp.F));
			bw.write("\t");
			bw.write(Double.toString(snp.p));
			bw.newLine();
			bw.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
}
