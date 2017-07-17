package net.maizegenetics.jGLiM.dm;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;

public class ModelEffectUtils {

	private ModelEffectUtils() {}
	
	public static DoubleMatrix getXtY(ModelEffect X, ModelEffect Y) {
		if (X instanceof FactorModelEffect) {
			FactorModelEffect fme = (FactorModelEffect) X;
			if (Y instanceof FactorModelEffect) {
				return fme.getXtX2((FactorModelEffect) Y);
			} else if (Y instanceof CovariateModelEffect) {
				return fme.getXty(((CovariateModelEffect) Y).getCovariate());
			} else if (Y instanceof NestedCovariateModelEffect) {
				return fme.getX().mult(Y.getX(), true, false);
			}
			
		} else if (X instanceof CovariateModelEffect) {
			double[] cov = ((CovariateModelEffect) X).getCovariate();
			if (Y instanceof FactorModelEffect) {
				return getXtY(Y,X).transpose();
			} else if (Y instanceof CovariateModelEffect) {
				return Y.getXty(cov);
			} else if (Y instanceof NestedCovariateModelEffect) {
				return Y.getXty(cov).transpose();
			}
			
		} else if (X instanceof NestedCovariateModelEffect) {
			if (Y instanceof FactorModelEffect) {
				return X.getX().mult(Y.getX(), true, false);
			} else if (Y instanceof CovariateModelEffect) {
				return X.getXty(((CovariateModelEffect) Y).getCovariate());
			} else if (Y instanceof NestedCovariateModelEffect) {
				return ((NestedCovariateModelEffect) X).getXtX2((NestedCovariateModelEffect) Y);
			}
			
		}
		return null;
	}
	
    public static int[] getIntegerLevels(Object[] originalLevels) {
    	int nLevels = originalLevels.length;
    	int[] intLevels = new int[nLevels];
    	HashMap<Object, Integer> levelMap = new HashMap<Object, Integer>();
    	for (int i = 0; i < nLevels; i++) {
    		Integer ndx = levelMap.get(originalLevels[i]);
    		if (ndx == null) {
    			ndx = new Integer(levelMap.size());
    			levelMap.put(originalLevels[i], ndx);
    		}
    		intLevels[i] = ndx.intValue();
    	}
    	return intLevels;
    }
    
    public static <T> int[] getIntegerLevels(T[] originalLevels, ArrayList<T> ids) {
    	int nLevels = originalLevels.length;
    	int[] intLevels = new int[nLevels];
    	HashMap<T, Integer> levelMap = new HashMap<T, Integer>();
    	for (int i = 0; i < nLevels; i++) {
    		Integer ndx = levelMap.get(originalLevels[i]);
    		if (ndx == null) {
    			ndx = new Integer(levelMap.size());
    			levelMap.put(originalLevels[i], ndx);
    		}
    		intLevels[i] = ndx.intValue();
    	}
    	
    	if (ids != null) {
        	TreeSet<Entry<T,Integer>> sortedEntries = new TreeSet<Entry<T,Integer>>(new Comparator<Entry<T,Integer>>(){

    			@Override
    			public int compare(Entry<T, Integer> arg0, Entry<T, Integer> arg1) {
    				return arg0.getValue().compareTo(arg1.getValue());
    			}
        		
        	});
        	
        	sortedEntries.addAll(levelMap.entrySet());
        	for (Entry<T, Integer> entry:sortedEntries) {
        		ids.add(entry.getKey());
        	}
    	}
    	
    	return intLevels;
    }
    
    public static <T> int[] getIntegerLevels(ArrayList<T> originalLevels, ArrayList<T> ids) {
    	int[] intLevels = new int[originalLevels.size()];
    	HashMap<T, Integer> levelMap = new HashMap<T,Integer>();
    	int count = 0;
    	for (T level:originalLevels) {
    		Integer ndx = levelMap.get(level);
    		if (ndx == null) {
    			ndx = new Integer(levelMap.size());
    			levelMap.put(level, ndx);
    		}
    		intLevels[count++] = ndx.intValue();
    	}

    	if (ids != null) {
        	TreeSet<Entry<T,Integer>> sortedEntries = new TreeSet<Entry<T,Integer>>(new Comparator<Entry<T,Integer>>(){

    			@Override
    			public int compare(Entry<T, Integer> arg0, Entry<T, Integer> arg1) {
    				return arg0.getValue().compareTo(arg1.getValue());
    			}
        		
        	});
        	
        	sortedEntries.addAll(levelMap.entrySet());
        	for (Entry<T, Integer> entry:sortedEntries) {
        		ids.add(entry.getKey());
        	}
    	}
    	return intLevels;
    }
    
    public static <T> int[] getIntegerLevels(ArrayList<T> originalLevels) {
    	return getIntegerLevels(originalLevels, null);
    }
    

}
