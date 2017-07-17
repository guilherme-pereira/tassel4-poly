package net.maizegenetics.gwas.modelfitter;

import java.util.ArrayList;
import java.util.LinkedList;

import net.maizegenetics.jGLiM.LinearModelUtils;
import net.maizegenetics.jGLiM.dm.CovariateModelEffect;
import net.maizegenetics.jGLiM.dm.FactorModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffectUtils;
import net.maizegenetics.jGLiM.dm.NestedCovariateModelEffect;
import net.maizegenetics.jGLiM.dm.PartitionedLinearModel;
import net.maizegenetics.jGLiM.dm.SweepFastLinearModel;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.MarkerPhenotypeAdapter;
import net.maizegenetics.pal.alignment.MarkerPhenotypeAdapterUtils;
import net.maizegenetics.pal.report.SimpleTableReport;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

public class StepwiseOLSModelFitter {
	private MarkerPhenotypeAdapter myData;
	
	//settable parameters
	private double[] enterlimits = null;
	private double[] exitlimits = null;
	private double enterlimit = 1e-5;
	private double exitlimit = 2e-5;
	private int maxNumberOfMarkers = 1000;
	private FactorModelEffect nestingEffect;
	private int nestingFactorIndex;
	private boolean isNested;
	private ArrayList<String> nestingFactorNames;

	//global variables used by the analysis
	private ArrayList<ModelEffect> currentModel;
	private int currentPhenotypeIndex;
	private int numberOfBaseModelEffects;
	private double[] y;
	private boolean[] missing;
	private ArrayList<String[]> factorList;
	private ArrayList<double[]> covariateList;
	private LinkedList<Object[]> resultRowsAnova = new LinkedList<Object[]>();
	private LinkedList<Object[]> rowsSiteEffectTable = new LinkedList<Object[]>();
    private MODEL_TYPE modelType = MODEL_TYPE.pvalue;
        
	private String datasetName;
        private double globalbestbic = Double.MAX_VALUE; //Rationale: the optimal model has minimum BIC. Thus, initialize bestbic with the largest possible double java can handle.
        private double globalbestmbic = Double.MAX_VALUE;
        public enum MODEL_TYPE {pvalue, bic, mbic};
        private boolean test;
        
        
	private final String[] anovaReportHeader = new String[]{"Trait", "Name","Locus","Position","df","SS","MS", "F", "pr>F"};
	
        public void setModelType(MODEL_TYPE modelType){
            this.modelType = modelType;
        }
        
	public StepwiseOLSModelFitter(MarkerPhenotypeAdapter anAdapter, String datasetName) {
		myData = anAdapter;
		this.datasetName = datasetName;
	}
	
	public DataSet runAnalysis() {
        //numberof markers
        int numberOfMarkers = myData.getNumberOfMarkers();
        
        //numbers of various model components
        int numberOfCovariates = myData.getNumberOfCovariates();
        int numberOfFactors = myData.getNumberOfFactors();
        int numberOfPhenotypes = myData.getNumberOfPhenotypes();
        
        //cycle through the phenotypes
        //notation: 
        //X is the design matrix without the markers, rows of X will be deleted if marker data is missing
        //Xm is the design matrix with markers
        //y is the data

        for (int ph = 0; ph < numberOfPhenotypes; ph++) {
        	currentPhenotypeIndex = ph;
        	if (enterlimits != null) enterlimit = enterlimits[ph];
        	if (exitlimits != null) exitlimit = exitlimits[ph];
        	
            //get phenotype data
            double[] phenotypeData = myData.getPhenotypeValues(ph);
            
            //keep track of missing rows
            missing = myData.getMissingPhenotypes(ph);

            //get factors
            factorList = MarkerPhenotypeAdapterUtils.getFactorList(myData, ph, missing);

            //get covariates
            covariateList = MarkerPhenotypeAdapterUtils.getCovariateList(myData, ph, missing);

            //remove missing values from the arrays
            int numberNotMissing = 0;
            int totalNumber = missing.length;
            for (boolean b:missing) if (!b) numberNotMissing++;
            
            int ptr = 0;
            y = new double[numberNotMissing];
            for (int i = 0; i < totalNumber; i++) {
            	if (!missing[i]) y[ptr++] = phenotypeData[i];
            }
            
            if (factorList != null) {
                int n = factorList.size();
                for (int f = 0 ; f < n; f++) {
                	String[] newfactor = new String[numberNotMissing];
                	String[] oldfactor = factorList.get(f);
                	ptr = 0;
                    for (int i = 0; i < totalNumber; i++) {
                    	if (!missing[i]) newfactor[ptr++] = oldfactor[i];
                    }
                    factorList.set(f, newfactor);
                }
            }
            
            if (covariateList != null) {
                int n = covariateList.size();
                for (int f = 0 ; f < n; f++) {
                	double[] newcov = new double[numberNotMissing];
                	double[] oldcov = covariateList.get(f);
                    for (int i = 0; i < totalNumber; i++) {
                    	if (!missing[i]) newcov[ptr++] = oldcov[i];
                    }
                    covariateList.set(f, newcov);
                }
            }
            fitModel();
        }

        return null;
	}
	
	public void fitModel() {
		//build the base model
		currentModel = new ArrayList<ModelEffect>();
		int numberOfTaxa = y.length;
		int[] mean = new int[numberOfTaxa];
		
		FactorModelEffect meanEffect = new FactorModelEffect(mean, false);
		meanEffect.setID("mean");
		currentModel.add(meanEffect);
		
		//add the factor effects
		if (factorList != null) {
			for (int f = 0; f < factorList.size(); f++) {
				ArrayList<String> ids = new ArrayList<String>();
				int[] levels = ModelEffectUtils.getIntegerLevels(factorList.get(f), ids);
				FactorModelEffect fme = new FactorModelEffect(levels, true, new Object[]{myData.getFactorName(f), ids});
				currentModel.add(fme);
				if (isNested && f == nestingFactorIndex) {
					nestingEffect = fme;
					nestingFactorNames = ids; 
				}
			}
		}
		
		//add the covariate effects
		if (covariateList != null) {
			for (int c = 0; c < covariateList.size(); c++) {
				CovariateModelEffect cme = new CovariateModelEffect(covariateList.get(c), myData.getCovariateName(c));
				currentModel.add(cme);
			}
		}
		numberOfBaseModelEffects = currentModel.size();
		
		while(forwardStep()) {
			while(backwardStep());
		}
		
		SweepFastLinearModel theLinearModel = new SweepFastLinearModel(currentModel, y);
		appendAnovaResults(theLinearModel);
		appendSiteEffectEstimates(theLinearModel);
		
	}
	
	public boolean forwardStep() {
		double bestss = 0;
                double bestbic = 1.7E+308;
                double bestmbic = 1.7E+308;
		ModelEffect besteffect = null;
		
		SweepFastLinearModel sflm = new SweepFastLinearModel(currentModel, y);
		PartitionedLinearModel plm = new PartitionedLinearModel(currentModel, sflm);
		int numberOfSites = myData.getNumberOfMarkers();
		
		for (int s = 0; s < numberOfSites; s++) {
			//create the appropriate marker effect
			ModelEffect markerEffect = null;
			SNP snp = new SNP(myData.getMarkerName(s), new Locus(myData.getLocusName(s)), (int) myData.getMarkerChromosomePosition(s), s);
			Object[] markerValues = myData.getMarkerValue(currentPhenotypeIndex, s);
			int n = markerValues.length;
                        
			if (myData.isMarkerDiscrete(s)) {
				//int n = markerValues.length;
				ArrayList<Object> markerIds = new ArrayList<Object>();
				int[] levels = ModelEffectUtils.getIntegerLevels(markerValues, markerIds);
				snp.alleles = markerIds;
				
				if (isNested) {
					//not implemented yet
					markerEffect = new FactorModelEffect(levels, true, snp);
				} else {
					markerEffect = new FactorModelEffect(levels, true, snp);
				}
			} else {
				//int n = markerValues.length;
				double[] cov = new double[n];
				for (int i = 0; i < n; i++) cov[i] = ((Double) markerValues[i]).doubleValue();
				
				if (isNested) {
					CovariateModelEffect cme = new CovariateModelEffect(cov);
					markerEffect = new NestedCovariateModelEffect(cme, nestingEffect);
					markerEffect.setID(snp);
				} else {
					markerEffect = new CovariateModelEffect(cov, snp);
				}
			}
			
			plm.testNewModelEffect(markerEffect);
			double modelss = plm.getModelSS();
                        
                        //Calculate the BIC
                        double errorss = plm.getErrorSS();
                        double modeldf = plm.getModeldf();
                        double p = modeldf+2; // modeldf corresponds to the slope estimates of each predictor, and 2 is for the intercept and the residual variance
                        double bic = (n*Math.log(errorss/n)) + (p*Math.log(n)) ; 
                                
                        //Calculate the mBIC
                        int numberOfTwoWayInteractions = (numberOfSites*(numberOfSites-1))/2;
                        double pForMbic = plm.getModeldf();
                        double qForMbic = 0; //This is going to be the number of two-way interaction terms. For now, this is not implemented, and thus I am setting it equal to zer
                        double mbic = (n*Math.log(errorss)) + ((pForMbic+qForMbic)*Math.log(n)) + (2*pForMbic*(Math.log((numberOfSites/2.2)-1))) + 
                                       (2*qForMbic*(Math.log((numberOfTwoWayInteractions/2.2)-1) ));
                        
                        switch(modelType){
                            case pvalue:
                                test = modelss > bestss;
                                break;
                            case bic:
                                test = bic < bestbic;
                                break;
                            case mbic:
                                //Not implemented yet
                                break;
                        }
			if (test) {
                                bestbic = bic;
				bestss = modelss;
				besteffect = markerEffect;
			}
		}
		
		//if the p-value for the select SNP is less than the enter limit, add it to the model and recalculate the model solution
		plm.testNewModelEffect(besteffect);
		double[] Fp = plm.getFp();
		if(modelType == MODEL_TYPE.pvalue){
                   if ( Fp[1] < enterlimit) {
			currentModel.add(besteffect);
			if (currentModel.size() == maxNumberOfMarkers + numberOfBaseModelEffects) return false;
			return true;
                    } else { 
			return false; 
                    }     
                } else if(modelType == MODEL_TYPE.bic){
                    if(bestbic < globalbestbic){
                        globalbestbic = bestbic;
                        currentModel.add(besteffect);
                        if (currentModel.size() == maxNumberOfMarkers + numberOfBaseModelEffects) return false;
                        return true;
                    } else{
                        return false;
                    } 
                } else if(modelType == MODEL_TYPE.mbic){
                    return false;
                }  else{
                    return false;
                }  
                   
	}
	
	public boolean backwardStep() {
		int numberOfTerms = currentModel.size();
		if (numberOfTerms <= numberOfBaseModelEffects) return false;
		
		SweepFastLinearModel sflm = new SweepFastLinearModel(currentModel, y);
	
		//find the model term (snps only) with the largest p-value
		double maxp = 0;
		double minF= -1;
		int maxterm = 0;
		double[] errorssdf = sflm.getResidualSSdf();

		for (int t = numberOfBaseModelEffects; t < numberOfTerms; t++) {
			double[] termssdf = sflm.getIncrementalSSdf(t);
			double F = termssdf[0]/termssdf[1]/errorssdf[0]*errorssdf[1];
			double p;
			try {
				p = LinearModelUtils.Ftest(F, termssdf[1], errorssdf[1]);
				if (p > maxp) {
					maxterm = t;
					maxp = p;
					minF = F; 
				}
			} catch(Exception e){
				p = Double.NaN;
			}
		}

		//if that maxp is >= exitlimit, then remove maxterm from the model, recalculate the model, and return true;
		if (maxp >= exitlimit) {
			ModelEffect me = currentModel.remove(maxterm);
			return true;
		}

		return false;
	}
	
	public LinkedList<Object[]> createReportRowsFromCurrentModel(SweepFastLinearModel sflm) {
		String traitname = myData.getPhenotypeName(currentPhenotypeIndex);
		int ncol = anovaReportHeader.length;
		LinkedList<Object[]> reportTable = new LinkedList<Object[]>();
		double[] residualSSdf = sflm.getResidualSSdf();
		int effectPtr = 0;
		for (ModelEffect me : currentModel) {
			Object[] reportRow = new Object[ncol];
			int ptr = 0;
			reportRow[ptr++] = traitname;
			reportRow[ptr++] = me.getID().toString();
			if (me.getID() instanceof SNP) {
				SNP snp = (SNP) me.getID();
				reportRow[ptr++] = snp.locus.getName();
				reportRow[ptr++] = Integer.toString(snp.position);
			} else {
				reportRow[ptr++] = "--";
				reportRow[ptr++] = "--";
			}
			double[] effectSSdf = sflm.getMarginalSSdf(effectPtr);
			double ms = effectSSdf[0] / effectSSdf[1];
			double Fval = ms / residualSSdf[0] * residualSSdf[1];
			double pval;
			try {
				pval = LinearModelUtils.Ftest(Fval, effectSSdf[1], residualSSdf[1]);
			} catch (Exception e) {
				pval = Double.NaN;
			}
			reportRow[ptr++] = new Integer((int) effectSSdf[1]);
			reportRow[ptr++] = new Double(effectSSdf[0]);
			reportRow[ptr++] = new Double(ms);
			reportRow[ptr++] = new Double(Fval);
			reportRow[ptr++] = new Double(pval);
			reportTable.add(reportRow);
			effectPtr++;
		}
		int ptr = 0;
		Object[] reportRow = new Object[ncol];
		reportRow[ptr++] = traitname;
		reportRow[ptr++] = "Error";
		reportRow[ptr++] = "--";
		reportRow[ptr++] = "--";
		reportRow[ptr++] = new Integer((int) residualSSdf[1]);
		reportRow[ptr++] = new Double(residualSSdf[0]);
		reportRow[ptr++] = new Double(residualSSdf[0]/residualSSdf[1]);
		reportRow[ptr++] = new Double(Double.NaN);
		reportRow[ptr++] = new Double(Double.NaN);
		reportTable.add(reportRow);
		
		return reportTable;
	}
	
	public Datum createReportFromCurrentModel(SweepFastLinearModel sflm) {
		String traitname = myData.getPhenotypeName(currentPhenotypeIndex);
		LinkedList<Object[]> reportTable = createReportRowsFromCurrentModel(sflm);
		Object[][] table = new Object[reportTable.size()][];
		reportTable.toArray(table);
		String reportName = "ANOVA for " + traitname + ", " + datasetName;
		TableReport tr = new SimpleTableReport(reportName, anovaReportHeader, table);
		return new Datum(reportName, tr, "");
	}
	
	public void appendAnovaResults(SweepFastLinearModel sflm) {
		resultRowsAnova.addAll(createReportRowsFromCurrentModel(sflm));
	}
	
	public TableReport getAnovaReport() {
		String reportName = "ANOVA table for " + datasetName;
		Object[][] table = new Object[resultRowsAnova.size()][];
		resultRowsAnova.toArray(table);
		return new SimpleTableReport(reportName, anovaReportHeader, table);
	}

	public TableReport getMarkerEffectReport() {
		if (rowsSiteEffectTable.size() == 0) return null;
		//columns are trait, snp name, locus, position, factor value, estimate
		String reportName = "Marker effects for " + datasetName;
		String[] reportHeader = new String[]{"Trait","Snp","Locus","Position","Within","Estimate"};
		Object[][] table = new Object[rowsSiteEffectTable.size()][];
		rowsSiteEffectTable.toArray(table);
		return new SimpleTableReport(reportName, reportHeader, table);
	}
	
	public void appendSiteEffectEstimates(SweepFastLinearModel sflm) {
		
		//columns are trait, snp name, locus, position, factor value, estimate
		int nBaseModelParameters = 0;
		String traitname = myData.getPhenotypeName(currentPhenotypeIndex);
		for (int i = 0; i < numberOfBaseModelEffects; i++) {
			nBaseModelParameters += currentModel.get(i).getEffectSize();
		}
		
		double[] beta = sflm.getBeta();
		int parameterIndex = nBaseModelParameters;
		for (int s = numberOfBaseModelEffects; s < currentModel.size(); s++) {
			ModelEffect me = currentModel.get(s);
			if (me instanceof NestedCovariateModelEffect) {
				NestedCovariateModelEffect ncme = (NestedCovariateModelEffect) me;
				FactorModelEffect fme = ncme.getFactorModelEffect();
				SNP snp = (SNP) ncme.getID();
				int n = nestingFactorNames.size();
				for (int i = 0; i < n; i++) {
					Object[] rowValues = new Object[6];
					int ptr = 0;
					rowValues[ptr++] = traitname;
					rowValues[ptr++] = snp.name;
					rowValues[ptr++] = snp.locus.getName();
					rowValues[ptr++] = new Integer(snp.position);
					rowValues[ptr++] = nestingFactorNames.get(i);
					rowValues[ptr++] = new Double(beta[parameterIndex++]);
					rowsSiteEffectTable.add(rowValues);
				}
			} else if (me instanceof CovariateModelEffect) {
				SNP snp = (SNP) me.getID();
				Object[] rowValues = new Object[6];
				int ptr = 0;
				rowValues[ptr++] = traitname;
				rowValues[ptr++] = snp.name;
				rowValues[ptr++] = snp.locus.getName();
				rowValues[ptr++] = new Integer(snp.position);
				rowValues[ptr++] = "NA";
				rowValues[ptr++] = new Double(beta[parameterIndex++]);
				rowsSiteEffectTable.add(rowValues);
			} else if (me instanceof FactorModelEffect) {
				
			}
		}
	}
	
	public void setEnterlimits(double[] enterlimits) {
		this.enterlimits = enterlimits;
	}

	public void setExitlimits(double[] exitlimits) {
		this.exitlimits = exitlimits;
	}

	public void setEnterlimit(double enterlimit) {
		this.enterlimit = enterlimit;
	}

	public void setExitlimit(double exitlimit) {
		this.exitlimit = exitlimit;
	}

	public void setMaxNumberOfMarkers(int maxNumberOfMarkers) {
		this.maxNumberOfMarkers = maxNumberOfMarkers;
	}

	public void setNestingEffectIndex(int nestingFactorIndex) {
		this.nestingFactorIndex = nestingFactorIndex;
	}

	public void setNested(boolean isNested) {
		this.isNested = isNested;
	}
}
