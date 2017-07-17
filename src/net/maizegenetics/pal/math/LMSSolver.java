// SUPGMABase.java
//
// (c) 1999-2003 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.math;

/**
 * <p>Title: LMSSolver </p>
 * <p>Description: An interface for objects that can perform a Least Mean Squares type regression.  </p>
 * @author Matthew Goode
 * @version 1.0
 */

public interface LMSSolver {

	public double[] solve(double[][] xMatrix, double[] dMatrix);

// -=-=-==-=-=-=-=-=-=-=--=-=-=-==--==-=--=-==-=-=-=-=-=--==-=--=
	/**
	 *
	 * <p>Title: Utils</p>
	 * A store of LMSSolver related utiltiies
	 */
	public static final class Utils {
		/**
		 * Obtain a simple LMSSolver that uses the "traditional" method for LMS stuff (IE, using row reduction to find the inverse,  with poor numerical accuracy)
		 * @return A LMSSolver object
		 */
		public static final LMSSolver getSimpleSolver() {
		  return SimpleSolver.INSTANCE;
		}
		// -==-=-=--==-
		private static final	class SimpleSolver implements LMSSolver {
			public static final LMSSolver INSTANCE = new SimpleSolver();
			public double[] solve( double[][] xMatrix, double[] dMatrix ) {
				Matrix m = new Matrix( xMatrix );
				Matrix mTranspose = m.getTranspose();
				Matrix ls = ( mTranspose.getMultiplied( m ) ).getInverse().getMultiplied( mTranspose );
				Matrix d = new Matrix( new double[][] {dMatrix} ).getTranspose();
				Matrix result = ls.getMultiplied( d );
				return result.toArray();
			}
		}
	}
}