package net.maizegenetics.baseplugins.genomicselection;
    
import net.maizegenetics.stats.EMMA.EMMAforDoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.*;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory.FactoryType;



public class RegRidgeEmmaDoubleMatrix {
		
	 
	private double[] GEBVs;
	private double[] mrkEsts;
	private DoubleMatrix GEBVsDM;
	private DoubleMatrix mrkEstsDM;
	
	private DoubleMatrix pheno;
	private boolean[] phenoMissing;
	private DoubleMatrix fixed;
	private DoubleMatrix geno;
	private int nLines;
	private int nObs;
	
	private DoubleMatrix kin = null;
	
	
	/** Constructor using native double arrays. Need to specify DoubleMatrix FactoryType
	 * @param phenotype; column matrix of phenotype values. Missing phenotype are coded as NaN. GEBVs are calculate for all lines including missing. 
	 * @param fixedEffects; 2D matrix of fixed effects; must not be singular
	 * @param genotypes; matrix of genotypes coded -1, 0, 1.  Missing Genotypes as NaN
	 * @param type; Tassel DoubleMatrix Type to use
	 */
	public RegRidgeEmmaDoubleMatrix(double[] phenotype, double[][] fixedEffects, double[][] genotypes, FactoryType type){
		DoubleMatrixFactory MF = new DoubleMatrixFactory(type);
		pheno = MF.make(phenotype.length, 1, phenotype);
		fixed = MF.make(fixedEffects);
		geno = MF.make(genotypes);
		phenoMissing = this.determineMissingPhenotypes();
	}
	
	
	/** Constructor using Tassel DoubleMatrix.
	 * @param phenotype; column matrix of phenotype values. Missing phenotype are coded as NaN. GEBVs are calculate for all lines including missing.
	 * @param fixedEffects; 2D matrix of fixed effects; must not be singular
	 * @param genotypes; matrix of genotypes coded -1, 0, 1.  Missing Genotypes as NaN
	 */
	public RegRidgeEmmaDoubleMatrix(DoubleMatrix phenotype, DoubleMatrix fixedEffects, DoubleMatrix genotypes){
		pheno = phenotype;
		fixed = fixedEffects;
		geno = genotypes;
		phenoMissing = this.determineMissingPhenotypes();
	}
	
	
	/** Constructor using Tassel DoubleMatrix. Use when kinship matrix from markers has already been calculated (XXt). 
	 * @param phenotype; column matrix of phenotype values. Missing phenotype are coded as NaN. GEBVs are calculate for all lines including missing.
	 * @param fixedEffects; 2D matrix of fixed effects; must not be singular
	 * @param genotypes; matrix of genotypes coded -1, 0, 1.  Missing Genotypes as NaN
	 * @param kinship; kinship matrix (XXt) calculated from marker matrix (X) 
	 */
	public RegRidgeEmmaDoubleMatrix(DoubleMatrix phenotype, DoubleMatrix fixedEffects, DoubleMatrix genotypes, DoubleMatrix kinship){
		pheno = phenotype;
		fixed = fixedEffects;
		geno = genotypes;
		kin = kinship;
		phenoMissing = this.determineMissingPhenotypes();
	}
	

	
	
	/**
	 * Solve the mixed model for line GEBVs and marker estimates
	 */
	public void solve() {
		
		double start = System.currentTimeMillis()/1000;
		
		int notMissing = 0;
		for (boolean b : phenoMissing) if (!b) notMissing++;
		int[] nmIndex = new int[notMissing];
		int count = 0;
		for (int i = 0; i < phenoMissing.length; i++) {
			if (!phenoMissing[i]) {nmIndex[count] = i; count++;}
		}
	
		DoubleMatrix nmPheno = pheno.getSelection(nmIndex, null);
		DoubleMatrix nmFixed = fixed.getSelection(nmIndex, null);

		nLines = pheno.numberOfRows();
		nObs = nmPheno.numberOfRows();
		
		this.replaceNaNwithMean();
		
		DoubleMatrix nmGeno = geno.getSelection(nmIndex, null);
		
		//calculate marker based kinship
		if (kin == null){kin = geno.mult(geno, false, true);}
		DoubleMatrix nmKin = kin.getSelection(nmIndex, nmIndex);
		
		// use EMMA to solve mixed model
		
		EMMAforDoubleMatrix lm = new EMMAforDoubleMatrix(nmPheno, nmFixed, nmKin, 0);
		lm.solve();
		double delta = lm.getDelta();
		System.out.println("Delta: " + delta);

	
		// solve mixed model for marker BLUPs
		DoubleMatrix aXZ = nmFixed.concatenate(nmGeno, false);
		 
		int nFixed = nmFixed.numberOfColumns();
		int nRandom = nmGeno.numberOfColumns();
		
		DoubleMatrix aCoefMat = aXZ.crossproduct();
		
		for (int i=nFixed; i<nFixed+nRandom; i++){
			aCoefMat.set(i,i,aCoefMat.get(i, i)+ delta); // TODO add weighting here. Currently weighted by including same line more than once
		}
		
		DoubleMatrix aXZty = aXZ.mult(nmPheno, true, false);
		
		int[] idx = new int[nRandom];// index for random effects
		for (int i=0; i<nRandom; i++){idx[i]=i+nFixed;}
		
		mrkEstsDM = aCoefMat.inverse().mult(aXZty, false, false).getSelection(idx, null);
		GEBVsDM = geno.mult(mrkEstsDM, false, false);
		
		mrkEsts = new double[nRandom];
		GEBVs = new double[nLines];
		
		for (int i=0; i<nRandom; i++){mrkEsts[i] = mrkEstsDM.get(i,0);}
		for (int i=0; i<nLines; i++){GEBVs[i] = GEBVsDM.get(i, 0);}
		
		double time = (System.currentTimeMillis()/1000 - start);
		System.out.println("Ran rrEMMA in " + time + " seconds");
		System.out.println();
	}


	
	/**
	 * @return double array of genomic estimated breeding values (GEBVs) for all lines
	 */
	public double[] getBlups(){
		return GEBVs;
	}
	
	/**
	 * @return double array of estimated marker effects for all markers
	 */
	public double[] getMrkBlups(){
		return mrkEsts;
	}
	
	/**
	 * @return Tassel DoubleMatrix of genomic estimated breeding values (GEBVs) for all lines
	 */
	public DoubleMatrix getGEBVsAsDoubleMatrix(){
		return GEBVsDM;
	}
	
	/**
	 * @return Tassel DoubleMatrix of estimated marker effects for all markers
	 */
	public DoubleMatrix getMrkEstsAsDoubleMatrix(){
		return mrkEstsDM;
	}
	
	
	/**
	 * replaces missing markers with the average marker score
	 */
	private void replaceNaNwithMean() {
	
		for (int c = 0; c < geno.numberOfColumns(); c++) {
			
		int nm = 0;
		double sum = 0;
		double mean;
		
			for (int r=0; r<nLines; r++){
				double val = geno.get(r, c);
				if (!Double.isNaN(val)) {
					nm++;
					sum += val;
				}
			}
			
			mean = sum / nm;

			for (int r = 0; r < nLines; r++) {
				if (Double.isNaN(geno.get(r, c))) { geno.set(r, c, mean); }
			}
		}

	}

	/**
	 * generates boolean array of missing phenotypes
	 */
	private boolean[] determineMissingPhenotypes(){
		boolean[] tmp = new boolean [pheno.numberOfRows()];
		for(int i = 0; i < pheno.numberOfRows(); i++ ) {
			if(Double.isNaN(pheno.get(i, 0))) {
				tmp[i] = true;
			} else {
				tmp[i] = false;
			}
		}
		return tmp;
	}
	
	
	
}
