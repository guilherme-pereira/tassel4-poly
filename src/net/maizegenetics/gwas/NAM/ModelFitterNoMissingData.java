package net.maizegenetics.gwas.NAM;

import java.io.BufferedWriter;
import java.util.ArrayList;

import net.maizegenetics.jGLiM.dm.FactorModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffectUtils;
import net.maizegenetics.jGLiM.dm.LinearModelForStepwiseRegression;
import net.maizegenetics.jGLiM.dm.ModelEffect;

public class ModelFitterNoMissingData extends ModelFitter {
	
	ArrayList<Integer> snpList = null;
	ArrayList<ModelEffect> effects;
	
	public ModelFitterNoMissingData(int chromosome, FileNames files) {
		super(chromosome, files);
	}
	
	public void testSnps() {

		//create the base model
		
		LinearModelForStepwiseRegression lmsr = getBaseModel();
		
		//setup the output
		BufferedWriter bw;
		bw = openOutputFile(files.chrsteps);
		
		StringBuilder sb = new StringBuilder("Position\tF\tp\tlog(1/p)");
		sb.append("\t").append("B97");
		sb.append("\t").append("CML103");
		sb.append("\t").append("CML228");
		sb.append("\t").append("CML247");
		sb.append("\t").append("CML277");
		sb.append("\t").append("CML322");
		sb.append("\t").append("CML333");
		sb.append("\t").append("CML52");
		sb.append("\t").append("CML69");
		sb.append("\t").append("Hp301");
		sb.append("\t").append("Il14H");
		sb.append("\t").append("Ki11");
		sb.append("\t").append("Ki3");
		sb.append("\t").append("Ky21");
		sb.append("\t").append("M162W");
		sb.append("\t").append("M37W");
		sb.append("\t").append("Mo17");
		sb.append("\t").append("Mo18W");
		sb.append("\t").append("MS71");
		sb.append("\t").append("NC350");
		sb.append("\t").append("NC358");
		sb.append("\t").append("Oh43");
		sb.append("\t").append("Oh7B");
		sb.append("\t").append("P39");
		sb.append("\t").append("Tx303");
		sb.append("\t").append("Tzi8");
		

		writeToOutput(sb.toString(), bw);
//		int snpcount = 0;
		while (snpdata.next()) {
//			if (snpcount++ % 1000 == 0) System.out.println("snpcount = " + snpcount);
			double[] parents = snpdata.getGenotype();
			int pos = snpdata.getPosition();
			
			//build and solve the model
			double[] snp = projectSnp(parents, pos, popIndex);
			double modelss = lmsr.testNewEffect(snp);
			double[] Fp = lmsr.getFpFromModelSS(modelss);
			
			sb = new StringBuilder();
			sb.append(pos);
			sb.append("\t").append(Fp[0]);
			sb.append("\t").append(Fp[1]);
			sb.append("\t").append(Math.log10(1/Fp[1]));
			
			for (int i = 0; i < 26; i++) {
				sb.append("\t");
				sb.append(parents[i]);
//				if (parents[i]) {
//					sb.append('0');
//				} else {
//					sb.append('1');
//				}
			}
			
			writeToOutput(sb.toString(), bw);
		}
		
		closeOutput(bw);
		
	}

	public LinearModelForStepwiseRegression getBaseModel() {
		int chromosome = files.chromosome;
		int nSamples = residualSamples.length;
		int[] mean = new int[nSamples];
		
		//set up the mean and pop model effects
		ArrayList<ModelEffect> effects = new ArrayList<ModelEffect>();
		ModelEffect memean = new FactorModelEffect(mean, false);
		effects.add(memean);
		ModelEffect mepop = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(pops), true);
		effects.add(mepop);
		
		//get data for this chromosome 
		double[] y = new double[nSamples];
		for (int i = 0; i < nSamples; i++) y[i] = residuals[i][chromosome - 1];
		
		//create the base model
		return new LinearModelForStepwiseRegression(effects, y);
	}

}
