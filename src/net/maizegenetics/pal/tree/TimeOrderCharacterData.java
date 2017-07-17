// TimeOrderCharacterData.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.tree;

import net.maizegenetics.pal.ids.IdGenerator;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.math.MersenneTwisterFast;
import net.maizegenetics.pal.util.BranchLimits;
import net.maizegenetics.pal.util.HeapSort;

import java.io.Serializable;

/**
 * Character data that expresses an order through time.
 *
 * @version $Id: TimeOrderCharacterData.java,v 1.3 2009/07/07 16:19:37 tcasstevens Exp $
 *
 * @author Alexei Drummond
 */
public class TimeOrderCharacterData implements Serializable, BranchLimits, UnitsProvider, IdGroup {

	//
	// Protected Stuff
	//

	/** Order of times */
	protected int[] timeOrdinals = null; //Is serialized

	/** Actual times of each sample */
	protected double[] times = null; //Is serialized

	/** the identifier group */
	protected IdGroup taxa; //Is serialized

	protected int units = Units.GENERATIONS; //Is serialized
	protected SubgroupHandler[] subgroups_;

	//
	// PRIVATE STUFF
	//

	/** Name of this character data */
	private String name = "Time/order character data";	//Is serialized

	//
	// Serialization code
	//
	private static final long serialVersionUID= 7672390862755080486L;

	//serialver -classpath ./classes net.maizegenetics.pal.tree.TimeOrderCharacterData
	private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
		out.writeByte(2); //Version number
		out.writeObject(timeOrdinals);
		out.writeObject(times);
		out.writeObject(taxa);
		out.writeInt(units);
		out.writeObject(name);
		out.writeObject(subgroups_);
	}

	private void readObject(java.io.ObjectInputStream in) throws java.io.IOException, ClassNotFoundException{
		byte version = in.readByte();
		switch(version) {
			case 1 : {
				timeOrdinals = (int[])in.readObject();
				times = (double[])in.readObject();
				taxa = (IdGroup)in.readObject();
				units = in.readInt();
				name = (String)in.readObject();
				break;
			}
			default : {
				timeOrdinals = (int[])in.readObject();
				times = (double[])in.readObject();
				taxa = (IdGroup)in.readObject();
				units = in.readInt();
				name = (String)in.readObject();
				subgroups_ = (SubgroupHandler[])in.readObject();
				break;
			}
		}
	}

	/**
	 * Parameterless constructor for superclasses.
	 */
	protected TimeOrderCharacterData() {}

	/**
		 * Clones a TimeOrderCharacterData object
		 * but the identifier positions match that of base (ie whichIdNumber(Name) returns the same as for base)
		 */
	private TimeOrderCharacterData(TimeOrderCharacterData toCopy, IdGroup base) {
		int size = toCopy.getIdCount();
		this.timeOrdinals = new int[size];
		final boolean hasTimes = toCopy.hasTimes();
		this.times = hasTimes ? new double[size] : null;
		for(int i = 0 ; i < size ; i++) {
			String name = toCopy.getIdentifier(i).getName();
			int baseLocation = base.whichIdNumber(name);
			if(baseLocation<0) {
				throw new IllegalArgumentException("Base does not contain:"+name);
			}
			this.timeOrdinals[baseLocation] = toCopy.timeOrdinals[i];
			if(hasTimes) {
				this.times[baseLocation] = toCopy.times[i];
			}
		}
		this.subgroups_ = SubgroupHandler.getCopy(toCopy.subgroups_,toCopy,base);
		this.units = toCopy.units;
		this.taxa = SimpleIdGroup.getInstance(base);
	}

	/**
	 * Constructor taking only IdGroup.
	 * Beware! This constructor does not initialize
	 * any time ordinals or times.
	 * @param taxa the taxa that this time data relates to.
	 * @param units the units of the times.
	 */
	public TimeOrderCharacterData(IdGroup taxa, int units) {
		this(taxa, units, false);
	}

	/**
	 * Constructor taking only IdGroup.
	 * @param taxa the taxa that this time data relates to.
	 * @param units the units of the times.
	 * @param contemp if true, all times are set to zero, else
	 * times are not set.
	 */
	public TimeOrderCharacterData(IdGroup taxa, int units, boolean contemp) {
		this.taxa = taxa;
		this.units = units;

		if (contemp) {
			double[] times = new double[taxa.getIdCount()];
			setTimes(times, units);
		}
	}

	/**
	 * Constructs a TimeOrderCharacterData with a number of
	 * equal-sized, evenly-spaced sampling times.
	 * @param numSeqsPerSample the number of taxa/sequences per sample time.
	 * @param numSamples the number of sample times.
	 * @param timeBetweenSamples the time between each pair of consecutive samples.
	 * @param units the units in which the times are expressed.
	 */
	public TimeOrderCharacterData(int numSeqsPerSample, int numSamples,
		double timeBetweenSamples, int units) {

		int n = numSeqsPerSample * numSamples;

		taxa = IdGenerator.createIdGroup(n);

		// create times and ordinals
		timeOrdinals = new int[taxa.getIdCount()];
		times = new double[taxa.getIdCount()];

		int index = 0;
		for (int i = 0; i < numSamples; i++) {
			for (int j = 0; j < numSeqsPerSample; j++) {
				times[index] = timeBetweenSamples * (double)i;
				timeOrdinals[index] = i;
				index += 1;
			}
		}

		this.units = units;
	}

	/**
	 * Returns a clone of the specified TimeOrderCharacterData
	 */
	public static TimeOrderCharacterData clone(TimeOrderCharacterData tocd) {
		return tocd.subset(tocd);
	}

	/**
	 * Extracts a subset of a TimeOrderCharacterData.
	 */
	public TimeOrderCharacterData subset(IdGroup staxa) {

		TimeOrderCharacterData subset =
			new TimeOrderCharacterData(staxa, getUnits());

		subset.timeOrdinals = new int[staxa.getIdCount()];
		if (hasTimes()) {
			subset.times = new double[staxa.getIdCount()];
		}

		for (int i = 0; i < subset.timeOrdinals.length; i++) {
			int index = taxa.whichIdNumber(staxa.getIdentifier(i).getName());
			subset.timeOrdinals[i] = timeOrdinals[index];

			if (hasTimes()) {
				subset.times[i] = times[index];
			}
		}
		return subset;
	}

	public int getUnits() {
		return units;
	}

	/**
	 * A means for define a subgroup.
	 * @param subgroups an array of integer arrays. Each array holds the indexes
	 * of the members that for that subgroup
	 */
	public final void setSubgroups(final int[][] subgroups) {
		this.subgroups_ = SubgroupHandler.create(subgroups);
	}
	/**
	 * A means for define a subgroup.
	 * @param subgroups an array indexes
	 * of the members that for the subgroup
	 */
	public final void setSubgroup(final int[] subgroup) {
		this.subgroups_ = SubgroupHandler.create(subgroup);
	}
	/**
	 * A means for define a subgroup.
	 * @param subgroups an array of names that represent the members of the subgroup.Non existent memebers are ignored
	 */
	public final void setSubgroup(final String[] subgroup) {
		this.subgroups_ = SubgroupHandler.create(this, subgroup);
	}
	/**
	 * A means for define subgroups.
	 * @param subgroups an array of String arrays. Each array holds the members
	 * for a particular subgroup. Nonexistent members are ignored.
	 * @note members can appear in more than one subgroup
	 */
	public final void setSubgroups(final String[][] subgroups) {
		this.subgroups_ = SubgroupHandler.create(this, subgroups);
	}

	public final boolean hasSubgroups() {
		return this.subgroups_!=null;
	}
	public final int getNumberOfSubgroups() {
		return (subgroups_==null? 0 : subgroups_.length);
	}
	/**
	 * Creates a TimeOrderCharacterData which is a subset of this sub group.
	 * Different subgroups may contain the same members
	 */
	public final TimeOrderCharacterData createSubgroup(int subgroupNumber) {
		return this.subgroups_[subgroupNumber].generateNewTOCD(this);
	}
	public final Identifier[] getSubgroupMembers(int subgroupNumber) {
		return this.subgroups_[subgroupNumber].getSubgroupMembers(this);
	}


	/**
	 * Sets the times, and works out what the ordinals should be.
	 */
	public void setTimes(double[] times, int units) {
		setTimes(times, units, true);
	}

	/**
	 * Sets the times.
	 * @param recalculateOrdinals true if ordinals should be
	 * recalculated from the times.
	 */
	public void setTimes(double[] times, int units, boolean recalculateOrdinals) {
		this.times = times;
		this.units = units;
		if (recalculateOrdinals) {
			setOrdinalsFromTimes();
		}
	}

	public TimeOrderCharacterData scale(double rate, int newUnits) {
		TimeOrderCharacterData scaled = clone(this);
		scaled.units = newUnits;
		for (int i = 0; i < times.length; i++) {
			scaled.times[i] = times[i] * rate;
		}
		return scaled;
	}
	/**
	 * Sets ordinals.
	 */
	public void setOrdinals(int[] ordinals) {
		timeOrdinals = ordinals;
	}
	/**
	 * @return the maximum time
	 */
	public double getMaximumTime() {
		if(times==null) {
			throw new RuntimeException("Error: getMaximumTime() called with no times");
		}
		double max = times[0];
		for(int i = 1 ; i < times.length ; i++) {
			if(times[i]>max) {
				max = times[i];
			}
		}
		return max;
	}
	/**
	 * @return the minimum time
	 */
	public double getMinimumTime() {
		if(times==null) {
			throw new RuntimeException("Error: getMinimumTime() called with no times");
		}
		double min = times[0];
		for(int i = 1 ; i < times.length ; i++) {
			if(times[i]<min) {
				min = times[i];
			}
		}
		return min;
	}

	/**
	 * Gets ordinals.
	 */
	public int[] getOrdinals() {
		return timeOrdinals;
	}

	/**
	 * Returns a copy of the times in the form of an array.
	 */
	public double[] getCopyOfTimes() {
		double[] copyTimes = new double[times.length];
		System.arraycopy(times, 0, copyTimes, 0, times.length);
		return copyTimes;
	}
	/**
	 * Creates a new TimeOrderCharacterData object with the same properites as this one
	 * but the identifier positions match that of base (ie whichIdNumber(Name) returns the same as for base)
	 * @throws IllegalArgumentException if the base ids don't match the ids of this tocd
	 *
	 */
	public TimeOrderCharacterData getReordered(IdGroup base) {
		return new TimeOrderCharacterData(this,base);
	}

	/**
	 * Remove time character data.
	 */
	public void removeTimes() {
		times = null;
	}
	/**
	 * Given an array of rates between samples (matching exactly the samples in order) then
	 * produces a TimeOrderCharacterData object that is timed by Expected Substitutions.
	 * Needs only sample information so no real time information required.
	 */
	public TimeOrderCharacterData generateExpectedSubsitutionsTimedTOCD(double[] sampleRates) {
		TimeOrderCharacterData result = new TimeOrderCharacterData(this,Units.EXPECTED_SUBSTITUTIONS);
		double[] times = new double[getIdCount()];
		for(int i = 0 ; i < times.length ; i++) {
			int sample = getTimeOrdinal(i);
			double total = 0;
			//Yes, I know this is inefficient but it's too late in the afternoon for me.
			for(int es = 0; es<sample ; es++) {	total+=sampleRates[es];	}
			times[i] = total;
		}
		result.setTimes(times,Units.EXPECTED_SUBSTITUTIONS);
		return result;
	}
	/**
	 * @return a dummy time order character data base on this that does have times, but
	 * the times are set to match the sample number. Eg. Ids in sample 0 have time 0, ids in sample 1 have time 1, and so on. The
	 * Units are UNKNOWN
	 */
	public TimeOrderCharacterData generateDummyTimedTOCD(double[] sampleRates) {
		TimeOrderCharacterData dummy = new TimeOrderCharacterData(this,Units.EXPECTED_SUBSTITUTIONS);
		double[] times = new double[getIdCount()];
		for(int i = 0 ; i < times.length ; i++) {
			times[i] = getTimeOrdinal(i);
		}
		dummy.setTimes(times, Units.UNKNOWN);
		return dummy;
	}



	/**
	 * Set time ordinals from another TimeOrderCharacterData.
	 * Select ordinals by matching names.
	 * @param tocd to take ordinals from.
	 */
	public void setOrdinals(TimeOrderCharacterData tocd) {
		setOrdinals(tocd, null, false);
	}

	public void setTimesAndOrdinals(TimeOrderCharacterData tocd) {
		setOrdinals(tocd, null, true);
	}

	/**
	 * Set time ordinals from another TimeOrderCharacterData.
	 * Select ordinals by matching names.
	 * @param tocd to take ordinals from
	 * @param idgroup use these labels to match indices in given tocd.
	 * @param doTimes if set then sets times as well
	 */
	public void setOrdinals(TimeOrderCharacterData tocd, IdGroup standard, boolean doTimes) {

		if (timeOrdinals == null) {
			timeOrdinals = new int[taxa.getIdCount()];
		}

		if (doTimes && tocd.hasTimes()) {
			times = new double[taxa.getIdCount()];
		}

		if (standard == null) {
			standard = tocd;
		}

		for (int i = 0; i < taxa.getIdCount(); i++) {

			String name = taxa.getIdentifier(i).getName();
			int index = standard.whichIdNumber(name);
			if (index == -1) {
				System.err.println("Identifiers don't match!");
				System.err.println("Trying to find: '" + name + "' in:");
				System.err.println(standard);
				//System.exit(1);
			}

			timeOrdinals[i] = tocd.getTimeOrdinal(index);
			if (doTimes && tocd.hasTimes()) {
				times[i] = tocd.getTime(index);
			}
		}
	}

	private final boolean equal(double a, double b, double epsilon) {
		return(Math.abs(a-b)<epsilon);
	}

	private void setOrdinalsFromTimes() {

		int[] indices = new int[times.length];
		timeOrdinals = new int[times.length];
		HeapSort.sort(times, indices);

		int ordinal = 0;
		int lastIndex = 0;
		timeOrdinals[indices[0]] = ordinal;

		for (int i = 1; i < indices.length; i++) {
			if (Math.abs(times[indices[i]] - times[indices[lastIndex]]) <= ABSTOL) {
				// this time is still within the tolerated error
			} else {
				// this is definitely a new time
				lastIndex = i;
				ordinal += 1;
			}
			timeOrdinals[indices[i]] = ordinal;
		}
	}

	/**
	 * Returns the number of characters per identifier
	 */
	public int getNumChars() {
		if (hasTimes()) {
			return 2;
		} else return 1;
	}

	/**
	 * Returns a name for this character data.
	 */
	public String getName() {
		return name;
	}

	/**
	 * Sets the name of this character data.
	 */
	public void setName(String name) {
		this.name = name;
	}

	public double getTime(int taxon) {
		return times[taxon];
	}
	/**
	 * Obtain the time of a particular ordinal
	 * @param ordinal The ordinal of interest
	 * @return The time of the input ordinal
	 * @throws IllegalArgumentException If no such ordinal
	 */
	public double getOrdinalTime(int ordinal) {
	  for(int i = 0 ; i < timeOrdinals.length ; i++) {
		  if(timeOrdinals[i]==ordinal) {
			  return times[i];
			}
		}
		throw new IllegalArgumentException("Unknown ordinal");
	}

	public double getTime(String taxonName) {
		int i = whichIdNumber(taxonName);
		return times[i];
	}


	/**
	 * NOTE: currently assumes times exist!
	 */
	public double getHeight(int taxon, double rate) {
		return times[taxon] * rate;
	}

	public int getTimeOrdinal(int taxon) {
		return timeOrdinals[taxon];
	}
	public int getTimeOrdinal(String taxonName) {
		int i = whichIdNumber(taxonName);
		return timeOrdinals[i];
	}
	public int getTimeOrdinal(Identifier taxonName) {
		int i = whichIdNumber(taxonName.getName());
		return timeOrdinals[i];
	}

	public boolean hasTimes() {
		return times != null;
	}

	/**
	 * Returns an ordered vector of unique times in this
	 * time order character data.
	 */
	public double[] getUniqueTimeArray() {
		int count = getSampleCount();
		double[] utimes = new double[count];
		for (int i = 0; i < times.length; i++) {
			utimes[getTimeOrdinal(i)] = times[i];
		}
		return utimes;
	}

	/**
	 * Returns a matrix of times between samples. A
	 * sample is any set of identifiers that have the same times.
	 */
	public double[][] getUniqueTimeMatrix() {
		double[] utimes = getUniqueTimeArray();
		int count = utimes.length;
		double[][] stimes = new double[count][count];
		for (int i = 0; i < count; i++) {
			for (int j = 0; j < count; j++) {
				stimes[i][j] = Math.abs(utimes[i] - utimes[j]);
			}
		}

		return stimes;
	}

	/**
	 * A sample is any set of identifiers that have the same times.
	 * @return the number of unique times in this data.
	 * @deprecated Use getOrdinalCount()
	 */
	public int getSampleCount() {
	  return getOrdinalCount();
	}
	/**
	 * @return the number of unique times in this data.
	 */
	public int getOrdinalCount() {
		int max = 0;
		for (int i = 0; i < timeOrdinals.length; i++) {
			if (timeOrdinals[i] > max) max = timeOrdinals[i];
		}
		return max + 1;
	}

	/**
	 * Returns a string representation of this time order character data.
	 */
	public String toString() {
		StringBuffer sb = new StringBuffer();

		sb.append("Identifier\t"+ (hasTimes() ? "Times\t" : "") + "Sample\n");
		for (int i = 0; i < taxa.getIdCount(); i++) {
			sb.append(taxa.getIdentifier(i) + "\t" +
				(hasTimes() ? getTime(i) + "\t" : "") +
				getTimeOrdinal(i)+"\n");
		}
		return new String(sb);
	}

	public void shuffleTimes() {
		MersenneTwisterFast mtf = new MersenneTwisterFast();

		int[] indices = mtf.shuffled(timeOrdinals.length);

		int[] newOrdinals = new int[timeOrdinals.length];
		double[] newTimes = null;
		if (hasTimes()) {
			newTimes = new double[times.length];
		}
		for (int i = 0; i < timeOrdinals.length; i++) {
			newOrdinals[i] = timeOrdinals[indices[i]];
			if (hasTimes()) { newTimes[i] = times[indices[i]]; }
		}

		timeOrdinals = newOrdinals;
		if (hasTimes()) times = newTimes;
	}

	//IdGroup interface
	public Identifier getIdentifier(int i) {return taxa.getIdentifier(i);}
	public void setIdentifier(int i, Identifier ident) { taxa.setIdentifier(i, ident); }
	public int getIdCount() { return taxa.getIdCount(); }
	public int whichIdNumber(String name) { return taxa.whichIdNumber(name); }
    public int whichIdNumber(Identifier id) { return taxa.whichIdNumber(id); }

	/**
	 * Return id group of this alignment.
	 * @deprecated TimeOrderCharacterData now implements IdGroup
	 */
	public IdGroup getIdGroup() { return taxa; }

	/**
	 * A simple utility method for generating a maximu mutation rate based
	 * on times. This is not guarranteed to be the best method. If the times are
	 * small (ie total differenct < 1), the maximum Mutation Rate is high. If the times are large
	 * (eg total differenct > 1) the maximum Mutation Rate is low. IE. the aim is to keep the Expected substitutions (mu*t) range reasonable (<1)
	 */
	public final double getSuggestedMaximumMutationRate() {
		double[] times = getUniqueTimeArray();
		double minDiff = times[1] - times[0];
		for(int i = 2 ; i < times.length ; i++) {
			double diff = times[i]-times[i-1];
			if(diff<minDiff) { minDiff = diff; }
		}
		return 5/minDiff;
	}

// ============================================================================
// ====== SubgroupHandler =====================================================
/**
 * Handles any subgroups
 */
	private static final class SubgroupHandler implements Serializable {
		private int[] subgroupIndexes_;

		//
		// Serialization code
		//
		private static final long serialVersionUID= 341384756632221L;

		//serialver -classpath ./classes net.maizegenetics.pal.tree.TimeOrderCharacterData
		private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
			out.writeByte(1); //Version number
			out.writeObject(subgroupIndexes_);
		}

		private void readObject(java.io.ObjectInputStream in) throws java.io.IOException, ClassNotFoundException{
			byte version = in.readByte();
			switch(version) {
				default : {
					subgroupIndexes_ = (int[])in.readObject();
					break;
				}
			}
		}
	// ===============
	// Constructors

		private SubgroupHandler(int[] subgroupIndexes) {
			this.subgroupIndexes_ = net.maizegenetics.pal.util.Utils.getCopy(subgroupIndexes);
		}
		private SubgroupHandler(SubgroupHandler toCopy, IdGroup oldBase, IdGroup newBase) {
			this(toCopy.subgroupIndexes_, oldBase,newBase);
		}
		private SubgroupHandler(int[] oldSubgroupIndexes, IdGroup oldBase, IdGroup newBase) {
			this.subgroupIndexes_ = new int[oldSubgroupIndexes.length];
			for(int i = 0 ; i < oldSubgroupIndexes.length ; i++) {
				String oldName = oldBase.getIdentifier(oldSubgroupIndexes[i]).getName();
				int newIndex = newBase.whichIdNumber(oldName);
				if(newIndex<0) {
					throw new IllegalArgumentException("Incompatible bases:"+oldName);
				}
				this.subgroupIndexes_[i] = newIndex;
			}
		}
		public Identifier[] getSubgroupMembers(TimeOrderCharacterData parent) {
			final int size = subgroupIndexes_.length;
			Identifier[] ids = new Identifier[size];
			for(int i = 0 ; i < size ; i++) {
				ids[i] = parent.getIdentifier(subgroupIndexes_[i]);
			}
			return ids;
		}
		public TimeOrderCharacterData generateNewTOCD(TimeOrderCharacterData parent) {
			final int size = subgroupIndexes_.length;
			Identifier[] ids = new Identifier[size];
			for(int i = 0 ; i < size ; i++) {
				ids[i] = parent.getIdentifier(subgroupIndexes_[i]);
			}
			TimeOrderCharacterData tocd = new TimeOrderCharacterData(new SimpleIdGroup(ids),parent.getUnits());

			if(parent.hasTimes()) {
				double[] times = new double[size];
				for(int i = 0 ; i < size ; i++) {
					times[i] = parent.getTime(subgroupIndexes_[i]);
				}
				tocd.setTimes(times,parent.getUnits());
			} else {
				int[] ordinals = new int[size];
				for(int i = 0 ; i < size ; i++) {
					ordinals[i] = parent.getTimeOrdinal(subgroupIndexes_[i]);
				}
				tocd.setOrdinals(ordinals);
			}
			return tocd;
		}
		public static final SubgroupHandler[] create(int[][] subgroupInfo) {
			SubgroupHandler[] handlers = new SubgroupHandler[subgroupInfo.length];
			for(int i = 0 ; i < handlers.length ; i++) {
				handlers[i] = new SubgroupHandler(subgroupInfo[i]);
			}
			return handlers;
		}
		public static final SubgroupHandler[] create(int[] subgroupInfo) {
			return new SubgroupHandler[] {
					new SubgroupHandler(subgroupInfo)
			};
		}
		private static final int[] getIndexes(final TimeOrderCharacterData parent, final String[] subgroupInfo) {
			int count = 0;
			for(int i = 0 ; i < subgroupInfo.length ; i++) {
				if(parent.whichIdNumber(subgroupInfo[i])>=0) {
					count++;
				}
			}
			final int[] indexes = new int[count];
			count = 0;
			for(int i = 0 ; i < subgroupInfo.length ; i++) {
				int parentIndex= parent.whichIdNumber(subgroupInfo[i]);
				if(parentIndex>=0) {
					indexes[count++] = parentIndex;
				}
			}
			System.out.println("Parent:"+parent);
			System.out.println("Indexes:"+net.maizegenetics.pal.util.Utils.toString(indexes));
			return indexes;
		}
		public static final SubgroupHandler[] create(TimeOrderCharacterData parent, String[][] subgroupInfo) {
			SubgroupHandler[] handlers = new SubgroupHandler[subgroupInfo.length];
			for(int i = 0 ; i < handlers.length ; i++) {
				handlers[i] = new SubgroupHandler(getIndexes(parent, subgroupInfo[i]));
			}
			return handlers;
		}
		/**
		 * Get a copy of the subgroup handlers such that the numbering is reorganised to match newbase
		 */
		public static final SubgroupHandler[] getCopy(SubgroupHandler[] toCopy, IdGroup oldBase, IdGroup newBase) {
			if(toCopy==null) {
				return null;
			}
			SubgroupHandler[] copy = new SubgroupHandler[toCopy.length];
			for(int i = 0 ; i < toCopy.length ; i++) {
				copy[i] = new SubgroupHandler(toCopy[i],oldBase,newBase);
			}
			return copy;
		}
		public static final SubgroupHandler[] create(TimeOrderCharacterData parent, String[] subgroupInfo) {
			return new SubgroupHandler[] {
					new SubgroupHandler(getIndexes(parent, subgroupInfo))
			};
		}
	}
}
