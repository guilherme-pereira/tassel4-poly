package net.maizegenetics.gwas.NAM;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;
import java.util.TreeMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import net.maizegenetics.jGLiM.LinearModelUtils;
import net.maizegenetics.jGLiM.dm.CovariateModelEffect;
import net.maizegenetics.jGLiM.dm.FactorModelEffect;
import net.maizegenetics.jGLiM.dm.LinearModelForStepwiseRegression;
import net.maizegenetics.jGLiM.dm.ModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffectUtils;
import net.maizegenetics.jGLiM.dm.SweepFastLinearModel;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;

public class ModelFitterBootstrapForward extends ModelFitter {
	static Random rng = new Random();
	ArrayList<int[]> populationArrays;
	int totalSubsampleSize;
	int[] popSubsampleSize;
	int[] subsample;
	double enterLimit = 1e-6;
	int maxThreads = 50;
	int numberOfThreads;
	int nSamples;
	
	public ModelFitterBootstrapForward(FileNames files) {
		super(files.chromosome, files);
		if (!Double.isNaN(files.enterlimit)) enterLimit = files.enterlimit;
	}

	@Override
	protected void testSnps() {
		long start = System.currentTimeMillis();
		System.out.println("Bootstrap analysis starting at " + System.currentTimeMillis());
		buildPopulationList();
		
		if (files.randomizeSnpOrder && files.threaded && files.fullModel) {
		       randomSnpThreadedAnalysis();
		} else {
		    if (files.threaded && files.fullModel) {
	            numberOfThreads = Runtime.getRuntime().availableProcessors();
	            threadedAnalysis();
		    } else if (files.fullModel) {
		    	simpleAnalysis();
		    } else {
		    	analyzeResiduals();
		    }
		}
		
		System.out.println("Bootstrap analysis finished, elapsed time = " + (System.currentTimeMillis() - start));
	}
	
	protected void simpleAnalysis() {
		//analyze bootstrap (subbagging) samples
		for (int i = 0; i < files.iterations; i++) {
			getSubsample();
			
			//fit forward regression
			LinearModelForStepwiseRegression lmsr = getBaseModel();
			
			SnpInfo nextSnp = findNextTerm(lmsr);
			int nsnpsAdded = 0;
			while (nextSnp.p < enterLimit && nsnpsAdded < files.maxsnps) {
				lmsr.addEffect(new CovariateModelEffect(sampleArray(nextSnp.genotype), nextSnp));
				nsnpsAdded++;
				nextSnp = findNextTerm(lmsr);
			}
			 
			try {
				recordResults(i + files.startIteration, lmsr);
			} catch (IOException e) {
				String msg = "Failed to write results for iteration " + i;
				e.printStackTrace();
			}
		}
	}
	
	protected void nonThreadedAnalysis() {
		int numberOfIterations = files.iterations;
		int chromosome = files.chromosome;
		ArrayList<ForwardRegressionModel> theModels = new ArrayList<ForwardRegressionModel>();
		for (int iter = 0; iter < numberOfIterations; iter++) {
			//get data for this chromosome 
			double[] data = new double[totalSubsampleSize];
			getSubsample();
			for (int i = 0; i < totalSubsampleSize; i++) data[i] = residuals[subsample[i]][chromosome - 1];
			
			theModels.add(new ForwardRegressionModelSubsample(getBaseEffects(), data, enterLimit, subsample, files.maxsnps));
		}

		boolean areAnyUpdateable = true;
		
		while (areAnyUpdateable) {
			snpdata.reset();
			boolean hasnext = snpdata.next();
			ArrayList<ForwardRegressionModel> updateableModels = new ArrayList<ForwardRegressionModel>();
			for (ForwardRegressionModel frm : theModels) {
				if (frm.wasUpdated) updateableModels.add(frm);
			}
			if (updateableModels.size() > 0) areAnyUpdateable = true;
			else areAnyUpdateable = false;
			
			while (hasnext) {
				// create a list of snps to test
				LinkedList<SnpInfo> snpList = new LinkedList<SnpInfo>();
				while (hasnext && snpList.size() <=500) {
					double[] snp = projectSnp(snpdata.getGenotype(), snpdata.getPosition(), popIndex);
					snpList.add(new SnpInfo(chromosome, snpdata.getPosition(),snpdata.getAllele(), snp, 1, 1));
					hasnext = snpdata.next();
				}

				//create and fire the updater threads
				ForwardRegressionUpdater theUpdater;
				theUpdater = new ForwardRegressionUpdater(updateableModels, snpList, !hasnext);
				theUpdater.run();
			}
		}
		
		//write the output
		int iterationCount = files.startIteration;
		for (ForwardRegressionModel frm : theModels) {
			
			try {
				recordResults(iterationCount++, frm.lmsr);
			} catch (IOException e) {
				System.err.println("Error recording results for iteration " + (iterationCount - 1));
				e.printStackTrace();
			}
		}
	}
	
	protected void threadedAnalysis() {
		int numberOfIterations = files.iterations;
		int chromosome = files.chromosome;
		ForwardRegressionModel[] theModels = new ForwardRegressionModel[numberOfIterations];
		ArrayList<ArrayList<ForwardRegressionModel>> theModelLists = new ArrayList<ArrayList<ForwardRegressionModel>>();
		for (int iter = 0; iter < numberOfIterations; iter++) {
			//get data for this chromosome 
			double[] data = new double[totalSubsampleSize];
			getSubsample();
			for (int i = 0; i < totalSubsampleSize; i++) data[i] = residuals[subsample[i]][chromosome - 1];
			
			theModels[iter] = new ForwardRegressionModelSubsample(getBaseEffects(), data, enterLimit, subsample, files.maxsnps);
		}


		boolean anyUpdateable = true;
		//create lists of models which need to have terms added
		//the number of lists equals the number of threads to be run
		//add snps then update the models; wash, rinse, repeat until all models are fit
		
		while (anyUpdateable) {
			System.out.println("Adding a term to the models.");
			long start = System.currentTimeMillis();
			theModelLists.clear();
			for (int i = 0; i < numberOfThreads; i++) {
				theModelLists.add(new ArrayList<ForwardRegressionModel>());
			}

			//create lists of forward regression models
			int currentList = 0;
			anyUpdateable = false;
			for (ForwardRegressionModel frm : theModels) {
				if (frm.wasUpdated) {
					theModelLists.get(currentList++).add(frm);
					if (currentList == numberOfThreads) currentList = 0;
					anyUpdateable = true;
				}
			}

			if (anyUpdateable) {
				snpdata.reset();
				boolean hasnext = snpdata.next();
				
				while (hasnext) {
					// create a list of snps to test
					LinkedList<SnpInfo> snpList = new LinkedList<SnpInfo>();
					while (hasnext && snpList.size() <=5000) {
						double[] snp = projectSnp(snpdata.getGenotype(), snpdata.getPosition(), popIndex);
						snpList.add(new SnpInfo(chromosome, snpdata.getPosition(),snpdata.getAllele(), snp, 1, 1));
						hasnext = snpdata.next();
					}
					
					//create and fire the updater threads
					ForwardRegressionUpdater[] updaters = null;
					updaters = new ForwardRegressionUpdater[numberOfThreads];
					for (int i = 0; i < numberOfThreads; i++) {
						updaters[i] = new ForwardRegressionUpdater(theModelLists.get(i), snpList, !hasnext);
						updaters[i].start();
					}
					
					//wait for the threads to finish
					for (int i = 0; i < numberOfThreads; i++) {
						try {
							updaters[i].join();
						} catch (InterruptedException e) {
						}
					}
				}
			}
			System.out.println("added a term in " + (System.currentTimeMillis() - start));
		}
		//write the output
		for (int i = 0; i < numberOfIterations; i++) {
			int thisIteration = files.startIteration + i;
			try {
				recordResults(thisIteration, theModels[i].lmsr);
			} catch (IOException e) {
				System.err.println("Error recording results for iteration " + thisIteration);
				e.printStackTrace();
			}
		}
	}
	
	protected void randomSnpThreadedAnalysis() {
		//Use a thread pool. For each thread:
		//Instantiate a forward regression model. Add that and a snpdata to a thread.
		//The thread should cycle through the snps several times until no more snps are added 
		//to the model. The finished model is a completed iteration, which is recorded.
		//The iterations can be numbered and written to output in the order finished.
		int chromosome = files.chromosome;
		//The thread pool
		int numberOfThreads = Runtime.getRuntime().availableProcessors();
		ExecutorService executor = Executors.newFixedThreadPool(numberOfThreads);
		LinkedList<Future<ForwardRegressionModel>> futureList = new LinkedList<Future<ForwardRegressionModel>>(); 
		
		//The loop for iterating through bootstraps
		for (int iter = 0; iter < files.iterations; iter++) {
			double[] data = new double[totalSubsampleSize];
			Integer[] pops = new Integer[totalSubsampleSize];
			int[] mean = new int[totalSubsampleSize];
			int[] sampleIndex = new int[totalSubsampleSize];
			getSubsample();
			for (int i = 0; i < totalSubsampleSize; i++) {
				data[i] = residuals[subsample[i]][chromosome - 1];
				pops[i] = new Integer(popIndex[subsample[i]]);
				sampleIndex[i] = subsample[i]; 
			}
			ArrayList<ModelEffect> initialEffects = new ArrayList<ModelEffect>();
			initialEffects.add(new FactorModelEffect(mean, false));
			initialEffects.add(new FactorModelEffect(ModelEffectUtils.getIntegerLevels(pops, null), true));
			ForwardRegressionModel frm = new ForwardRegressionModelSubsample(initialEffects, data, enterLimit, sampleIndex, files.maxsnps);
			futureList.add(executor.submit(new ForwardRegressionFitter(this, frm, snpdata.getCopy())));
		}
		System.out.println("All models have been created. Fitting will proceed.");
		executor.shutdown();
		
		//Output the results to a file
		int modelCount = files.startIteration;
		while (futureList.size() > 0) {
			Future<ForwardRegressionModel> f = futureList.remove();
			try {
				ForwardRegressionModel frm = f.get();
				recordResults(modelCount, frm.lmsr);
				
				System.out.println("Results recorded for iteration " + modelCount);
			} catch (InterruptedException e) {
				System.err.println("Iteration " + modelCount + " failed.");
				e.printStackTrace();
			} catch (ExecutionException e) {
				System.err.println("Iteration " + modelCount + " failed.");
				e.printStackTrace();
			} catch (IOException e) {
				System.err.println("Error writing results for iteration " + modelCount);
				e.printStackTrace();
			}
			modelCount++;
		}
	}

	protected void analyzeResiduals() {
		//analyze bootstrap (subbagging) samples
		System.out.println("analyze by residuals");
		
		//get data for this chromosome
		int chromosome = files.chromosome;
		
		for (int i = 0; i < files.iterations; i++) {
			getSubsample();
			double[] data = new double[totalSubsampleSize];
			for (int j = 0; j < totalSubsampleSize; j++) data[j] = residuals[subsample[j]][chromosome - 1];

			//fit forward regression
			LinearModelForStepwiseRegression lmsr = getBaseModel();
			SnpInfo nextSnp = findNextTerm(lmsr, true);
			int nsnpsAdded = 0;
			while (nextSnp.p < enterLimit && nsnpsAdded < files.maxsnps) {
				System.out.println("Adding " + nextSnp.pos + ", " + nextSnp.F + ", " + nextSnp.p);
				lmsr.addEffect(new CovariateModelEffect(sampleArray(nextSnp.genotype), nextSnp));
				nsnpsAdded++;
				if (nsnpsAdded < files.maxsnps) nextSnp = findNextTerm(lmsr, true);
			}
			System.out.println("-------------------------------------------------------------------------------");
			
			//sent results to output
			try {
				recordResults(i + files.startIteration, lmsr);
			} catch (IOException e) {
				String msg = "Failed to write results for iteration " + i;
				e.printStackTrace();
			}
		}		
	}
	
	//for the next term (max modelss) need pos, allele, projection, F, p
	protected SnpInfo findNextTerm(LinearModelForStepwiseRegression lmsr) {
		int chromosome = files.chromosome;
		SnpInfo bestSnp = new SnpInfo(chromosome, -1, "n/a", new double[]{}, 0,1);
		double bestSS = 0;
		snpdata.reset();
		while (snpdata.next()) {
			double[] parents = snpdata.getGenotype();
			int pos = snpdata.getPosition();
			
			//build and solve the model
			double[] snp = projectSnp(parents, pos, popIndex);
			
			double ms =  lmsr.testNewEffect(sampleArray(snp));
			
			if (ms > bestSS) {
				bestSS = ms;
				double[] Fp = lmsr.getFpFromModelSS(ms);
				bestSnp = new SnpInfo(chromosome, pos, snpdata.getAllele(), snp, Fp[0], Fp[1]);
			}
		}
		System.out.println(bestSnp.pos + ", " + bestSnp.allele + ", " + bestSnp.F + ", " + bestSnp.p);
		return bestSnp;
	}

	protected SnpInfo findNextTerm(LinearModelForStepwiseRegression lmsr, boolean useResidual) {
		DoubleMatrix residual = lmsr.getLinearModel().getResiduals();
		int chromosome = files.chromosome;
		double bestSS = 0;
		int bestPos = -1;
		String bestSnpAllele = "";
		double[] bestsnp = null;
		
		snpdata.reset();
		while (snpdata.next()) {
			double[] parents = snpdata.getGenotype();
			int pos = snpdata.getPosition();
			
			//build and solve the model
			double[] snp = projectSnp(parents, pos, popIndex);
			double SS = SnpSSUsingResidual(residual, sampleArray(snp));
			
			if (SS > bestSS) {
				bestSS = SS;
				bestPos = snpdata.getPosition();
				bestSnpAllele = snpdata.getAllele();
				bestsnp = snp;
			}
		}
		
		double ms =  lmsr.testNewEffect(sampleArray(bestsnp));
		double[] Fp = lmsr.getFpFromModelSS(ms);
		SnpInfo bestSnp = new SnpInfo(chromosome, bestPos, bestSnpAllele, bestsnp, Fp[0], Fp[1]);
//		System.out.println(bestSnp.pos + ", " + bestSnp.allele + ", " + bestSnp.F + ", " + bestSnp.p);
		return bestSnp;
	}
	
	public static synchronized void shuffle(int[] array) {
	    // i is the number of items remaining to be shuffled.
		int n = array.length;
	    for (int i = n; i > 1; i--) {
	        // Pick a random element to swap with the i-th element.
	        int j = rng.nextInt(i);  // 0 <= j <= i-1 (0-based array)
	        // Swap array elements.
	        int tmp = array[j];
	        array[j] = array[i-1];
	        array[i-1] = tmp;
	    }
	}

	private void buildPopulationList() {
		TreeMap<Integer, LinkedList<Integer>> popmap = new TreeMap<Integer, LinkedList<Integer>>();
		int count = 0;
		for (int pop : popIndex) {
			LinkedList<Integer> poplist = popmap.get(pop);
			if (poplist == null) {
				poplist = new LinkedList<Integer>();
				popmap.put(pop, poplist);
			}
			poplist.add(count++);
		}

		populationArrays = new ArrayList<int[]>();
		int n = popmap.size();
		totalSubsampleSize = 0;
		popSubsampleSize = new int[n];

		Iterator<LinkedList<Integer>> mit = popmap.values().iterator();
		int pcount = 0;
		while (mit.hasNext()) {
			LinkedList<Integer> plist = mit.next();
			int m = plist.size();
			int[] samples = new int[m];
			count = 0;
			for (Integer sample : plist) {
				samples[count++] = sample.intValue();
			}
			populationArrays.add(samples);

			popSubsampleSize[pcount] = m * 8 / 10;
			totalSubsampleSize += popSubsampleSize[pcount];
			pcount++;

		}

	}
	
	private void getSubsample() {
		int start = 0;
		int npops = populationArrays.size();
		subsample = new int[totalSubsampleSize];
		
		for (int i = 0; i < npops; i++) {
			int[] samples = populationArrays.get(i);
			shuffle(samples);
			System.arraycopy(samples, 0, subsample, start, popSubsampleSize[i]);
			start += popSubsampleSize[i];
		}
	}

	private double[] sampleArray(double[] in) {
		int nIn = in.length;
		int nOut = totalSubsampleSize;
		double[] out = new double[nOut];
		for (int i = 0; i < nOut; i++) {
			out[i] = in[subsample[i]];
		}
		return out;
	}
	
	protected LinearModelForStepwiseRegression getBaseModel() {
		ArrayList<ModelEffect> effects = getBaseEffects();
		int chromosome = files.chromosome;
		//get data for this chromosome 
		double[] data = new double[totalSubsampleSize];
		for (int i = 0; i < totalSubsampleSize; i++) data[i] = residuals[subsample[i]][chromosome - 1];
		
		//create the base model
		return new LinearModelForStepwiseRegression(effects, data);
	}

	protected ArrayList<ModelEffect> getBaseEffects() {
		int[] mean = new int[totalSubsampleSize];

		//set up the mean and pop model effects
		ArrayList<ModelEffect> effects = new ArrayList<ModelEffect>();
		ModelEffect memean = new FactorModelEffect(mean, false);
		effects.add(memean);
		ArrayList<Integer> populations = new ArrayList<Integer>();
		for (int sample : subsample) populations.add(popIndex[sample]);
		ModelEffect mepop = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(populations), true);
		effects.add(mepop);
		
		return effects;
	}
	
	protected double SnpSSUsingResidual(DoubleMatrix residual, double[] snp) {
		int n = residual.numberOfRows();
		double sumx = 0;
		double sumy = 0;
		double sumxsq = 0;
		double sumxy = 0;
		
		for (int i = 0; i < n; i++) {
			double x = snp[i];
			double y = residual.get(i, 0);
			sumx += x;
			sumy += y;
			sumxsq += x * x;
			sumxy += x * y;
		}
		
		double N = n;
		double snpSS = sumxsq - sumx / N * sumx;
		double sumprod = sumxy - sumx / N * sumy;
		
		double b = sumprod / snpSS;
		double yhatSS = b * sumprod;
		return yhatSS;
	}
	
	protected double SnpSSUsingResidual2(DoubleMatrix residual, double[] snp) {
		int n = residual.numberOfRows();
		DoubleMatrix snpMatrix = DoubleMatrixFactory.DEFAULT.make(n, 1, snp);
		double sumx = snpMatrix.columnSum(0);
		double sumy = residual.columnSum(0);
		double sumxsq = snpMatrix.crossproduct().get(0, 0);
		double sumxy = snpMatrix.crossproduct(residual).get(0, 0);
		
		double N = n;
		double snpSS = sumxsq - sumx / N * sumx;
		double sumprod = sumxy - sumx / N * sumy;
		
		double b = sumprod / snpSS;
		double yhatSS = b * sumprod;
		return yhatSS;
	}
	
	private void initializeOutput() throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(files.chrmodel));
		bw.write("chromosome\tposition\tcM\tallele\teffect\tF\tpvalue\titeration");
		bw.newLine();
		bw.close();
	}
	
	private void recordResults(int iteration, LinearModelForStepwiseRegression lmsr) throws IOException {
		int chromosome = files.chromosome;
		String tab = "\t";
		ArrayList<ModelEffect> effects = lmsr.getModelEffects();
		int nsnps = effects.size() - 2;
		double[] beta = lmsr.getLinearModel().getBeta();
		int nbeta = beta.length;
		int start = nbeta - nsnps;
		
		//write to model file
		BufferedWriter bw = new BufferedWriter(new FileWriter(files.chrmodel, true));
		double[] errorssdf = lmsr.getLinearModel().getResidualSSdf();
		for (int i = 0; i < nsnps; i++) {
			ModelEffect thiseffect = effects.get(i + 2);
			SnpInfo snpinfo = (SnpInfo) thiseffect.getID();
			
			//calculate F and p
			double[] snpssdf = lmsr.getLinearModel().getMarginalSSdf(i + 2);
			double F = snpssdf[0] / snpssdf[1] / errorssdf[0] * errorssdf[1];
			double p;
			try {
				p = LinearModelUtils.Ftest(F, snpssdf[1], errorssdf[1]);
			} catch (Exception e) {
				p = Double.NaN;
			}
			
			StringBuilder sb = new StringBuilder();
			sb.append(chromosome);
			sb.append(tab).append(snpinfo.pos);
			sb.append(tab).append(theAGPMap.getCmFromPosition(chromosome, snpinfo.pos));
			sb.append(tab).append(snpinfo.allele);
			sb.append(tab).append(beta[start + i]);
			sb.append(tab).append(F);
			sb.append(tab).append(p);
			sb.append(tab).append(iteration);
			bw.write(sb.toString());
			bw.newLine();
		}
		bw.close();
		
		//write to step file
		bw = new BufferedWriter(new FileWriter(files.chrsteps, true));
        for (int i = 0; i < nsnps; i++) {
            ModelEffect thiseffect = effects.get(i + 2);
            SnpInfo snp = (SnpInfo) thiseffect.getID();
            bw.write(Integer.toString(snp.chromosome));
            bw.write("\t");
            bw.write(Integer.toString(snp.pos));
            bw.write("\t");
            bw.write(snp.allele);
            bw.write("\t");
            bw.write(Double.toString(snp.F));
            bw.write("\t");
            bw.write(Double.toString(snp.p));
            bw.write("\t");
            bw.write(Integer.toString(iteration));
            bw.newLine();
        }
        bw.close();
		
	}

	private void recordResultsFromAnalyzeByResiduals(int iteration, LinearModelForStepwiseRegression lmsr) throws IOException {
		int chromosome = files.chromosome;
		String tab = "\t";
		ArrayList<ModelEffect> effects = lmsr.getModelEffects();
		int nsnps = effects.size() - 2;
		double[] beta = lmsr.getLinearModel().getBeta();
		int nbeta = beta.length;
		int start = nbeta - nsnps;
		
		//write to model file
		BufferedWriter bw = new BufferedWriter(new FileWriter(files.chrmodel, true));
		double[] errorssdf = lmsr.getLinearModel().getResidualSSdf();
		for (int i = 0; i < nsnps; i++) {
			ModelEffect thiseffect = effects.get(i + 2);
			SnpInfo snpinfo = (SnpInfo) thiseffect.getID();
			
			//calculate F and p
			double[] snpssdf = lmsr.getLinearModel().getMarginalSSdf(i + 2);
			double F = snpssdf[0] / snpssdf[1] / errorssdf[0] * errorssdf[1];
			double p;
			try {
				p = LinearModelUtils.Ftest(F, snpssdf[1], errorssdf[1]);
			} catch (Exception e) {
				p = Double.NaN;
			}
			
			StringBuilder sb = new StringBuilder();
			sb.append(chromosome);
			sb.append(tab).append(snpinfo.pos);
			sb.append(tab).append(theAGPMap.getCmFromPosition(chromosome, snpinfo.pos));
			sb.append(tab).append(snpinfo.allele);
			sb.append(tab).append(beta[start + i]);
			sb.append(tab).append(F);
			sb.append(tab).append(p);
			sb.append(tab).append(iteration);
			bw.write(sb.toString());
			bw.newLine();
		}
		bw.close();
		
		//write to step file
		bw = new BufferedWriter(new FileWriter(files.chrsteps, true));
        for (int i = 0; i < nsnps; i++) {
            ModelEffect thiseffect = effects.get(i + 2);
            SnpInfo snp = (SnpInfo) thiseffect.getID();
            bw.write(Integer.toString(snp.chromosome));
            bw.write("\t");
            bw.write(Integer.toString(snp.pos));
            bw.write("\t");
            bw.write(snp.allele);
            bw.write("\t");
            bw.write(Double.toString(snp.F));
            bw.write("\t");
            bw.write(Double.toString(snp.p));
            bw.write("\t");
            bw.write(Integer.toString(iteration));
            bw.newLine();
        }
        bw.close();
		
	}
}
