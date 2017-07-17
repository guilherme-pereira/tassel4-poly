// MutationRateModel.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.tree;

import net.maizegenetics.pal.io.FormattedOutput;
import net.maizegenetics.pal.math.OrthogonalHints;
import net.maizegenetics.pal.report.Report;
import net.maizegenetics.pal.util.Parameterized;

import java.io.Serializable;


/**
 * This abstract class contains methods that are of general use for
 * modelling mutation rate changes over time.
 *
 * @version $Id: MutationRateModel.java,v 1.1 2007/01/12 03:26:17 tcasstevens Exp $
 *
 * @author Alexei Drummond
 */
public abstract class MutationRateModel implements Units,
        Parameterized, Report, Cloneable, Serializable
{


	//
	// Private and protected stuff
	//

	protected FormattedOutput fo;

	/**
	 * Units in which time units are measured.
	 */
	private int units;

	private double maximumMutationRate_;


	private static final long serialVersionUID=-1755051453782951214L;

	//serialver -classpath ./classes net.maizegenetics.pal.mep.ConstantMutationRate

	private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
		out.writeByte(1); //Version number
		out.writeInt(units);
		out.writeDouble(maximumMutationRate_);
	}

	private void readObject(java.io.ObjectInputStream in) throws java.io.IOException, ClassNotFoundException{
		byte version = in.readByte();
		switch(version) {
			default : {
				units = in.readInt();
				maximumMutationRate_ = in.readDouble();
				fo = FormattedOutput.getInstance();
				break;
			}
		}
	}


	protected MutationRateModel(int units, double maximumMutationRate)
	{
		setUnits(units,maximumMutationRate);
		fo = FormattedOutput.getInstance();
	}
	protected MutationRateModel(MutationRateModel toCopy) {
		this.units = toCopy.units;
		this.maximumMutationRate_ = toCopy.maximumMutationRate_;
		fo = FormattedOutput.getInstance();
	}


	public abstract Object clone();

	//
	// Because I don't like casting if I don't have to
	//
	public abstract MutationRateModel getCopy();

	//
	// functions that define a mutation rate model (left for subclass)
	//

	/**
	 * Gets the mutation rate, value of mu(t) at time t.
	 */
	public abstract double getMutationRate(double t);

	/**
	 * Returns integral of mutation rate function
	 * (= integral mu(x) dx from 0 to t).
	 */
	public abstract double getExpectedSubstitutions(double t);

	/**
	 * Return the time at which expected substitutions has occurred.
	 */
	public double getTime(double expectedSubs) {
		return getEndTime(expectedSubs,0);
	}

	/**
	 * Return the end time at which expected substitutions has occurred, given we start at start time
	 */
	public abstract double getEndTime(double expectedSubs, double startTime);

	/**
	 * Linearly scales this mutation rate model.
	 * @param scale getExpectedSubstitutions should return scale instead of 1.0 at time t.
	 */
	public abstract void scale(double scale);

	// Parameterized and Report interface is also left for subclass


	// general functions

	/**
	 * Calculates the integral 1/mu(x) dx between start and finish.
	 */
	public double getExpectedSubstitutions(double start, double finish)
	{
		return getExpectedSubstitutions(finish) - getExpectedSubstitutions(start);
	}
	/**
	 * @throws IllegalArgumentException if units of this Model doenot match
	 * the units of the TimeOrderCharacterData object (toScale).
	 * @return a TimeOrderCharacterData scaled to use EXPECTED_SUBSTITUTIONS based on
	 * this MutationRateModel
	 */
	public TimeOrderCharacterData scale(TimeOrderCharacterData toScale) {
		if(getUnits()!=toScale.getUnits()) {
			throw new IllegalArgumentException("Incompatible units, expecting "+getUnits()+", found (in toScale) "+toScale.getUnits());
		}
		TimeOrderCharacterData scaled = toScale.clone(toScale);
		double[] times = new double[scaled.getIdCount()];
		for (int i = 0; i < times.length; i++) {
			times[i] = getExpectedSubstitutions(scaled.getTime(i));
		}
		scaled.setTimes(times,Units.EXPECTED_SUBSTITUTIONS,false);
		return scaled;
	}

	/**
	 * sets units of measurement.
	 * @throws IllegalArgumentException if units are ExpectedSubstitutions
	 *
	 * @param u units
	 * @param the maximumMutationRate that is allowable, given the units. This needs to be given intelligently.
	 */
	public final void setUnits(int u, double maximumMutationRate)
	{
		if(u==Units.EXPECTED_SUBSTITUTIONS) { throw new IllegalArgumentException("Units cannot be Expected Substitutions!"); }
		units = u;
		this.maximumMutationRate_ = maximumMutationRate;
	}
	/**
	 * @return the maximum mutation rate as indicated by the user
	 */
	protected final double getMaximumMutationRate() { return maximumMutationRate_; }


	/**
	 * returns units of measurement.
	 */
	public int getUnits()
	{
		return units;
	}

	/**
	 * Overide if there is any orthogonal hint information available
	 * @return null
	 */
	public OrthogonalHints getOrthogonalHints() {
		return null;
	}

	public abstract String toSingleLine();

	public abstract Factory generateFactory();
// ===========================================================================
// ==== Factory interface
	/**
	 * An interface for objects which generate fresh MutationRAteModels
	 */
	public interface Factory {
		/**
		 * Request a new MutationRateModel instance
		 */
		public MutationRateModel generateNewModel();
	}

}
