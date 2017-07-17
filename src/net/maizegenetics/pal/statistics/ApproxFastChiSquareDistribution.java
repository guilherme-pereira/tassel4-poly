// ChiSquareDistribution.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.statistics;

import java.util.Arrays;




/**
 * chi-square distribution
 * (distribution of sum of squares of n uniform random variables)
 *
 * (Parameter: n; mean: n; variance: 2*n)
 *
 * The chi-square distribution is a special case of the Gamma distribution
 * (shape parameter = n/2.0, scale = 2.0).
 *
 * @version $Id: ChiSquareDistribution.java,v 1.1 2007/01/12 03:26:16 tcasstevens Exp $
 *
 * @author Korbinian Strimmer
 */
public class ApproxFastChiSquareDistribution extends GammaDistribution
{
    double[][] precomputed;
    int maxN, maxX;

    public ApproxFastChiSquareDistribution(int maxX, int maxN) {
        this.maxN=maxN;
        this.maxX=maxX;
        precomputed=new double[maxN][maxX+1];
        Arrays.fill(precomputed[0], 0);
        for (int n = 1; n < maxN; n++) {
            for (int x = 0; x <= maxX; x++) {
                precomputed[n][x] = cdf(x,n);
   //             System.out.printf("%d %d %g %g %n",x, n, precomputed[n][x], cdfFastApprox(x,n));
            }  
        }
//        for (int n = 1; n < maxN; n+=20) {
//            for (double x = 0; x < maxX; x++) {
//               // precomputed[n][x] = cdf(x,n);
//                System.out.printf("%g %d %g %g %g %g %n",x+0.7, n, cdf(x+0.7,n), cdfFastApprox(x,n),
//                        cdfFastApprox(x+0.7,n), cdfFastApprox(x+1,n));
//            }  
//        }
    }
	//
	// Public stuff
	//
    

	/**
	 * probability density function of the chi-square distribution
	 * 
	 * @param x argument
	 * @param n degrees of freedom
	 *
	 * @return pdf value
	 */
	public static double pdf(double x, double n)
	{
		return pdf(x, n/2.0, 2.0);
	}

	/**
	 * cumulative density function of the chi-square distribution
	 * 
	 * @param x argument
	 * @param n degrees of freedom
	 *
	 * @return cdf value
	 */
	public static double cdf(double x, double n)
	{
		return cdf(x, n/2.0, 2.0);
	}
        
        /**
	 * cumulative density function of the chi-square distribution
	 * 
	 * @param x argument
	 * @param n degrees of freedom
	 *
	 * @return cdf value
	 */
	public double cdfFastApprox(double x, double n)
	{       int lx=(int)x;
                int ln=(int)n;
                int lxFloor=(int)Math.floor(x);
                int lxCeil=(int)Math.ceil(x);
                double d=x-(double)lxFloor;
                double p=((d*precomputed[ln][lxCeil])+((1.0-d)*precomputed[ln][lxFloor]));
                return p;
	//	return precomputed[ln][lx];
	}


	/**
	 * quantile (inverse cumulative density function) of the chi-square distribution
	 * 
	 * @param x argument
	 * @param n degrees of freedom
	 *
	 * @return icdf value
	 */
	public static double quantile(double y, double n)
	{
		return quantile(y, n/2.0, 2.0);
	}
	
	/**
	 * mean of the chi-square distribution
	 * 
	 * @param n degrees of freedom
	 *
	 * @return mean
	 */
	public static double mean(double n)
	{
		return n;
	}

	/**
	 * variance of the chi-square distribution
	 * 
	 * @param n degrees of freedom
	 *
	 * @return variance
	 */
	public static double variance(double n)
	{
		return 2.0*n;
	}
}
