package net.maizegenetics.gwas.NAM;

import java.util.List;

public class ForwardRegressionUpdater extends Thread {
	List<ForwardRegressionModel> myModels;
	List<SnpInfo> mySnps;
	boolean update;
	boolean areAnyUpdateable;
	
	/**
	 * This class creates a thread to test a series of SNPs in a series of models or to update the models.
	 * @param models	a List of ForwardRegressionModels objects
	 * @param snps	a List of SnpInfo objects. This can be null if update is true.
	 * @param update	if true, the models will be updated. Otherwise the snps will be tested in each model.
	 */
	public ForwardRegressionUpdater(List<ForwardRegressionModel> models, List<SnpInfo> snps, boolean update) {
		myModels = models;
		mySnps = snps;
		this.update = update;
	}
	
	@Override
	public void run() {
		for (SnpInfo snp : mySnps) {
			for (ForwardRegressionModel frm : myModels) {
				frm.addNextSnp(snp);
			}
		}
		
		if (update) {
			for (ForwardRegressionModel frm : myModels) {
					frm.updateModel();
			}
		} 
	}
	
	public boolean areAnyUpadateable() {return areAnyUpdateable;}
}
