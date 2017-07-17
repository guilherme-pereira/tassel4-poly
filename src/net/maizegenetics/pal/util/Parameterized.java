// Parameterized.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.util;




/**
 * interface for class with (optimizable) parameters
 *
 * @version $Id: Parameterized.java,v 1.1 2007/01/12 03:26:18 tcasstevens Exp $
 *
 * @author Korbinian Strimmer
 */
public interface Parameterized
{
	/**
	 * get number of parameters
	 *
	 * @return number of parameters
	 */
	int getNumParameters();

	/**
	 * set model parameter
	 *
	 * @param param  parameter value
	 * @param n  parameter number
	 */
	void setParameter(double param, int n);

	/**
	 * get model parameter
	 *
	 * @param n  parameter number
	 *
	 * @return parameter value
	 */
	double getParameter(int n);


	/**
	 * set standard errors for model parameter
	 *
	 * @param paramSE  standard error of parameter value
	 * @param n parameter number
	 */
	void setParameterSE(double paramSE, int n);


	/**
	 * get lower parameter limit
	 *
	 * @param n parameter number
	 *
	 * @return lower bound
	 */
	double getLowerLimit(int n);

	/**
	 * get upper parameter limit
	 *
	 * @param n parameter number
	 *
	 * @return upper bound
	 */
	double getUpperLimit(int n);


	/**
	 * get default value of parameter
	 *
	 * @param n parameter number
	 *
	 * @return default value
	 */
	double getDefaultValue(int n);

// ============================================================================

	/**
	 * A Utility class for using as the superclass to subclasses which work by adding functionality to
	 * a general Parameterized object (the base Parameterized object)
	 */
	abstract public static class ParameterizedUser //extends PalObjectListener.EventGenerator
     {
		private Parameterized base_;

		protected ParameterizedUser(Parameterized base) {
			this.base_ = base;
		}

		/**
		 * Sub class should call setParameterizedBase() at some point if using this constructor... otherwise many explosions and the flashing of lights.
		 */
		protected ParameterizedUser() {}

		/**
		 * Set's the base parameterized object...
		 */
		protected void setParameterizedBase(Parameterized base) {
			this.base_ = base;
		}

		public int getNumParameters() { return base_.getNumParameters(); }

		public void setParameter(double param, int n) {
			base_.setParameter(param,n);
			//fireParametersChangedEvent();
		}
		protected double getRandomParameterValue(int parameter) {
			double min = base_.getLowerLimit(parameter);
			double max = base_.getUpperLimit(parameter);
			return Math.random()*(max-min)+min;
		}
		protected void randomiseParameters() {
			final int numberParameters = base_.getNumParameters();
			for(int i = 0 ; i < numberParameters ; i++) {
				base_.setParameter(getRandomParameterValue(i),i);
			}
		}
		public double getParameter(int n) { return base_.getParameter(n); }

		public void setParameterSE(double paramSE, int n) { base_.setParameterSE(paramSE,n);  }

		public double getLowerLimit(int n) {	return base_.getLowerLimit(n);	}

		public double getUpperLimit(int n) { return base_.getUpperLimit(n); }
		public double getDefaultValue(int n) { return base_.getDefaultValue(n); }
	}


	/**
	 * A Utility class for using as the superclass to subclasses which are based on double arrays
	 */
	abstract public static class ParameterizedBase //extends PalObjectListener.EventGenerator
 {
		private double[] parameters_;
		private double[] lowerLimits_;
		private double[] upperLimits_;
		private double[] defaultValues_;
		private double[] parametersSE_;
		/**
		 * Builds a parameterized objects around a set of double arrays
		 * @note uses the actual passed objects, so subclasses can maintain references
		 * to the arrays and their values will be refected in the parameterized interface
		 */
		protected ParameterizedBase(
			double[] parameters,
			double[] lowerLimits,
			double[] upperLimits,
			double[] defaultValues
			) {
			setSource(parameters,lowerLimits,upperLimits,defaultValues);
		}
		/**
		 * Builds a parameterized objects around a set of double arrays
		 * @note uses the actual passed objects, so subclasses can maintain references
		 * to the arrays and their values will be refected in the parameterized interface
		 */
		protected ParameterizedBase(
			double[] parameters,
			double[] lowerLimits,
			double[] upperLimits,
			double[] defaultValues,
			double[] parametersSE
			) {
			setSource(parameters,lowerLimits,upperLimits,defaultValues, parametersSE);
		}
		/**
		 * Builds a parameterized objects around a set of double arrays. The user needs
		 * to call setSource() at some point to correctly set things up
		 * @note uses the actual passed objects, so subclasses can maintain references
		 * to the arrays and their values will be refected in the parameterized interface
		 */
		protected ParameterizedBase() { }
		/**
		 * Set's the base arrays...
		 * @throws IllegalArgumentException if arrays are not the same size
		 */
		protected void setSource(
			double[] parameters,
			double[] lowerLimits,
			double[] upperLimits,
			double[] defaultValues
		) {
			setSource(parameters,lowerLimits,upperLimits,defaultValues,new double[parameters.length]);
		}
		/**
		 * Set's the base arrays...
		 * @throws IllegalArgumentException if arrays are not the same size
		 */
		protected void setSource(
			double[] parameters,
			double[] lowerLimits,
			double[] upperLimits,
			double[] defaultValues,
			double[] parametersSE
		) {
			this.parameters_ = parameters;
			this.lowerLimits_ = lowerLimits;
			this.upperLimits_ = upperLimits;
			this.defaultValues_ = defaultValues;
			this.parametersSE_ = parametersSE;
			int length = parameters_.length;
			if(lowerLimits.length!=length) { sizeError();	}
			if(upperLimits.length!=length) { sizeError();	}
			if(defaultValues.length!=length) { sizeError();	}
			if(parametersSE.length!=length) { sizeError();	}
		}
		private final void sizeError() {
			throw new IllegalArgumentException("All arrays do not match in size");
		}

		public int getNumParameters() { return parameters_.length; }

		public void setParameter(double param, int n) {
			this.parameters_[n] = param;
			//fireParametersChangedEvent();
		}

		public double getParameter(int n) { return this.parameters_[n]; }

		public void setParameterSE(double paramSE, int n) { this.parametersSE_[n] = paramSE;  }

		public double getLowerLimit(int n) {	return lowerLimits_[n];	}
		public double getUpperLimit(int n) { return upperLimits_[n]; }
		public double getDefaultValue(int n) { return defaultValues_[n]; }
		protected double[] getParametersSE() { return parametersSE_; }
	}

	/**
	 * NullParameterized Object
	 * Can be used by subclasses to implement parameterized without actually having paramters
	 */
	abstract public static class Null implements Parameterized {
		protected Null() { }
		public int getNumParameters() { return 0; }
		private double error() { throw new RuntimeException("Assertion error : should not be called as no parameters!"); }
		public void setParameter(double param, int n) { error(); }
		public double getParameter(int n) { return error(); }
		public void setParameterSE(double paramSE, int n) { error(); }

		public double getLowerLimit(int n) {	return error();	}
		public double getUpperLimit(int n) { return error(); }
		public double getDefaultValue(int n) { return error(); }
	}

// ============================================================================
// ============ Utils ==========================================================
	public static final class Utils {
		/**
		 * @return the current paramters of source in a double array
		 */
		public static final double[] getParameters(Parameterized source) {
			double[] params = new double[source.getNumParameters()];
			for(int i = 0 ; i < params.length ; i++) {
				params[i] = source.getParameter(i);
			}
			return params;
		}
		/**
		 * Create a wrapper around a set of double arrays to create a parameterized object
		 * (changes to parameterized object change given arrays)
		 */
		public static final Parameterized createParametizedWrapper(double[] parameters, double[] lowerLimits, double[] upperLimits, double[] defaultValues) {
			return new ParameterizedWrapper(parameters,lowerLimits,upperLimits,defaultValues);
		}
		/**
		 * Create a wrapper around a set of double arrays to create a parameterized object
		 * (changes to parameterized object change given arrays)
		 */
		public static final Parameterized createParametizedWrapper(double[] parameters, double[] lowerLimits, double[] upperLimits, double[] defaultValues, double[] parametersSE) {
			return new ParameterizedWrapper(parameters,lowerLimits,upperLimits,defaultValues, parametersSE);
		}

		/**
		 * Combine multiple parameterized objects into one
		 */
		public static final Parameterized combine(Parameterized[] bases) {
			return new MultiParameterized(bases);
		}
		/**
		 * Combine multiple parameterized objects into one
		 */
		public static final Parameterized combine(Parameterized baseOne, Parameterized baseTwo) {
			return new MultiParameterized(new Parameterized[] {baseOne, baseTwo});
		}

		// ============================================================================
		// ============ ParameterizedWrapper ============================================
		private static final class ParameterizedWrapper extends ParameterizedBase implements Parameterized {
			public ParameterizedWrapper(double[] parameters, double[] lowerLimits, double[] upperLimits, double[] defaultValues) {
				super(parameters,lowerLimits,upperLimits,defaultValues);
			}
			public ParameterizedWrapper(double[] parameters,	double[] lowerLimits, double[] upperLimits,	double[] defaultValues,	double[] parametersSE) {
				super(parameters,lowerLimits,upperLimits,defaultValues, parametersSE);
			}

		}
	}

}
