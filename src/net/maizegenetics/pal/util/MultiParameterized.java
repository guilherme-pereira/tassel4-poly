// MultiParameterized.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.util;

/**
 * <p>Title: MultiParameterized</p>
 * <p>Description: A utility class for integrating multiple parameterized objects into one</p>
 * @author Matthew Goode
 * @version 1.0
 */

public class MultiParameterized implements Parameterized {
	private final Parameterized[] bases_;
	private final int[] baseLookup_;
	private final int[] baseParameter_;
	private final ParameterAccessWatcher watcher_;

	public MultiParameterized(Parameterized[] bases) {
		this(bases,null);
	}
	public MultiParameterized(Parameterized[] bases, ParameterAccessWatcher watcher) {
		this.watcher_ = watcher;
		final int numberOfParameters = countParameters(bases);

		this.baseLookup_ = new int[numberOfParameters];
		this.baseParameter_ = new int[numberOfParameters];
		this.bases_ = bases;
		setup(bases_, baseLookup_, baseParameter_);
	}
	protected static final int countParameters(final Parameterized[] ps) {
		int count = 0;
		for(int i = 0 ; i < ps.length ; i++) {
			count+=ps[i].getNumParameters();
		}
		return count;
	}
	protected static final void setup(final Parameterized[] ps, final int[] baseLookup, final int[] baseParameter) {
		int index = 0;
		for(int pIndex = 0 ; pIndex < ps.length ; pIndex++) {
			Parameterized p = ps[pIndex];
			final int numberOfParameters = p.getNumParameters();
			for(int i = 0 ; i < numberOfParameters ; i++) {
				baseLookup[index] = pIndex;
				baseParameter[index] = i;
				index++;
			}
		}
	}
	public final int getNumParameters() { return baseLookup_.length; }
	public void setParameterSE(final double paramSE, final int n) {	bases_[baseLookup_[n]].setParameterSE(paramSE, baseParameter_[n]);			}
	public void setParameter(final double param, final int n) {
		bases_[baseLookup_[n]].setParameter(param, baseParameter_[n]);
		if(watcher_!=null) {
			watcher_.parameterSet(bases_[baseLookup_[n]],param, baseParameter_[n]);
		}
	}

	public double getParameter(final int n) {		return bases_[baseLookup_[n]].getParameter(baseParameter_[n]);	}
	public double getUpperLimit(final int n) {	return bases_[baseLookup_[n]].getUpperLimit(baseParameter_[n]); 	}
	public double getLowerLimit(final int n) {	return bases_[baseLookup_[n]].getLowerLimit(baseParameter_[n]);			}
	public double getDefaultValue(final int n) {	return bases_[baseLookup_[n]].getDefaultValue(baseParameter_[n]);		}

	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append("Multi\n");
		for(int i = 0 ; i < bases_.length ;i++) {
			sb.append(i+":"+bases_[i]+"\n");
		}
		return sb.toString();
	}
	// -==--=-=-==--=-==-=-=-=--=
	/**
	 * An interface for classes that wish to find out about particular parameter access
	 */
	public static interface ParameterAccessWatcher {
		public void parameterSet(Parameterized baseParameterized, double param, int localParameter);
	}
}
