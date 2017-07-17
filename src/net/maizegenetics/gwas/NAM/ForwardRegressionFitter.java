package net.maizegenetics.gwas.NAM;

import java.util.concurrent.Callable;

public class ForwardRegressionFitter implements Callable<ForwardRegressionModel> {

	ForwardRegressionModel myModel;
	SnpData mySnpData;
	ModelFitter parentFitter;
	
	public ForwardRegressionFitter(ModelFitter parent, ForwardRegressionModel frm, SnpData snps) {
		myModel = frm;
		mySnpData = snps;
		parentFitter = parent;
		
	}

	@Override
	public ForwardRegressionModel call() {
		boolean isUpdateable = true;
		while (isUpdateable) {
			mySnpData.reset();

			//debug
//			System.out.println("SnpData has pointer set to " + ((SnpDataRandomOrder)mySnpData).currentIndex);
			while (mySnpData.next()) {
				//project the snp here
				int pos = mySnpData.getPosition();
				double[] snpScores = parentFitter.projectSnp(mySnpData.getGenotype(), pos);
				SnpInfo snp = new SnpInfo(mySnpData.chromosome, pos, mySnpData.getAllele(), snpScores, 1, 1);
				myModel.addNextSnp(snp);
			}
			isUpdateable = myModel.updateModel();
		}
		return myModel;
	}
	
}
