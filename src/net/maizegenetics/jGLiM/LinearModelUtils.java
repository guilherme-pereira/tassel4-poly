package net.maizegenetics.jGLiM;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import net.maizegenetics.jGLiM.dm.FactorModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffectUtils;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;

/**
 * @author Peter
 *
 */
public class LinearModelUtils {
	
	//prevents this class from being instantiated
	private LinearModelUtils(){
		
	}
	
    /**
     * Calculates the p-value associated with an F statistic.  The returned p-value is the probability
     * that a greater F is drawn from the F-distribution.
     *
     * @param F             - value of the F statistic
     * @param numeratordf   - degreees of freedom in the numerator of F
     * @param denominatordf - degreees of freedom in the denominator of F
     * @return p-value
     */
    public static double Ftest(double F, double numeratordf, double denominatordf) {
        double k = denominatordf / (denominatordf + numeratordf * F);
        return cern.jet.stat.Gamma.incompleteBeta(denominatordf / 2, numeratordf / 2, k);
    }

    /**
     * @param factorList	an ArrayList of String[], where each String[] contains the names of the levels of a factor
     * @param covariateList
     * @param missing
     * @return
     */
    public static DoubleMatrix createFixedEffectsArray(ArrayList<String[]> factorList, ArrayList<double[]> covariateList, boolean[] missing) {
    	int numberOfFactors;
    	int numberOfCovariates;
    	if (factorList == null) numberOfFactors = 0;
    	else numberOfFactors = factorList.size();
    	if (covariateList == null) numberOfCovariates = 0;
    	else numberOfCovariates = covariateList.size();
    	int numberOfEffects = 1 + numberOfFactors + numberOfCovariates;
    	DoubleMatrix[][] theMatrices = new DoubleMatrix[1][numberOfEffects];
    	int numberOfObs = 0;
    	for (boolean m:missing) if (!m) numberOfObs++;
    	
    	//the mean
    	int count = 0;
    	theMatrices[0][count++] = DoubleMatrixFactory.DEFAULT.make(numberOfObs, 1, 1.0);

    	for (int i = 0; i < numberOfFactors; i++) {
    		String[] nonMissingFactorLevels = getNonMissingElements(factorList.get(i), missing);
    		int[] levels = ModelEffectUtils.getIntegerLevels(nonMissingFactorLevels, null);
    		FactorModelEffect fme = new FactorModelEffect(levels, true);
    		theMatrices[0][count++] = fme.getX();
    	}

    	for (int i = 0; i < numberOfCovariates; i++) {
    		double[] nonMissingValues = getNonMissingElements(covariateList.get(i), missing);
    		theMatrices[0][count++] = DoubleMatrixFactory.DEFAULT.make(numberOfObs, 1, nonMissingValues);
    	}
    	
    	if (theMatrices[0].length == 1) return theMatrices[0][0];
    	return DoubleMatrixFactory.DEFAULT.compose(theMatrices);
    }
    
    /**
     * @param <T> 
     * @param array	an array of type T 
     * @param missing	an array of booleans equal to true if that element of the array should be deleted, false otherwise
     * @return	all the non-missing elements of array in the original order
     */
    public static <T> T[] getNonMissingElements(T[] array, boolean[] missing) {
    	int numberNotMissing = 0;
    	for (boolean m:missing) if (!m) numberNotMissing++;
    	T[] reducedArray = Arrays.copyOf(array, numberNotMissing);
    	int n = array.length;
    	int count = 0;
    	for (int i = 0; i < n; i++) {
    		if (!missing[i]) reducedArray[count++] = array[i];
    	}
    	return reducedArray;
    }

    /**
     * @param array	an array of doubles 
     * @param missing	an array of booleans equal to true if that element of the array should be deleted, false otherwise
     * @return	all the non-missing elements of array in the original order
     */
    public static double[] getNonMissingElements(double[] array, boolean[] missing) {
    	int numberNotMissing = 0;
    	for (boolean m:missing) if (!m) numberNotMissing++;
    	double[] reducedArray = new double[numberNotMissing++];
    	int n = array.length;
    	int count = 0;
    	for (int i = 0; i < n; i++) {
    		if (!missing[i]) reducedArray[count++] = array[i];
    	}
    	return reducedArray;
    }

	public static void shuffle(double[] source, Random randomizer) {
		int n = source.length;
		//the following algorithm from http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle (3/11/2011)
//		To shuffle an array a of n elements:
//			  for i from n - 1 down to 1 do
//			       j = random integer with 0 <= j <= i
//			       exchange a[j] and a[i]
		
		for (int i = n - 1; i > 0; i--) {
			int j = randomizer.nextInt(i + 1);
			double k = source[j];
			source[j] = source[i];
			source[i] = k;
		}
	}

}
