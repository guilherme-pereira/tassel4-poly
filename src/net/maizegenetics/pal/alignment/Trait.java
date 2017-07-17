package net.maizegenetics.pal.alignment;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

/**
 * Data descriptor that contains information about a trait, factor, covariate, 
 * or other item for which data may be stored. It may contain information 
 * about the experimental unit on which data was collected (factors). It may contain information 
 * about how the data is to be used in an analysis. Data values will be stored as doubles. 
 * In the case of discrete variables, the value stored will be 0, 1, ... indicating the variable level. 
 * For discrete variables, labels or names corresponding to each variable level would be stored as well. 
 * @author Peter Bradbury
 *
 */
/**
 * @author Peter
 *
 */
/**
 * @author Peter
 *
 */
public class Trait implements Comparable<Trait>, Serializable {
	protected String name;
	protected boolean isDiscrete;
	protected String type;
	protected HashMap<String, Object> propertyMap;

	public static final String TYPE_DATA = "data";
	public static final String TYPE_COVARIATE = "covariate";
	public static final String TYPE_FACTOR = "factor";
	public static final String TYPE_MARKER = "marker";
	public static final String TYPE_EXCLUDE = "exclude";
	
	public static final String PROP_LEVEL_LABELS = "level_labels";
	public static final String PROP_FACTORS = "factors";
	public static final String PROP_CHROMOSOME = "chromosome";
	public static final String PROP_POSITION = "position";
	public static final String PROP_LOCUS = "locus";
	public static final String PROP_LOCUS_POSITION = "locusposition";
	
	public static final String FACTOR_ENV = "env";

	/**
	 * Constructor
	 * @param name			this trait
	 * @param isDiscrete	true if this trait is discrete, false if it is continuous
	 * @param type			the trait type; data, covariate, factor, marker, or exclude
	 * @param properties	a property map, which can be null
	 */
	public Trait(String name, boolean isDiscrete, String type, Map<String, Object> properties) {
		this.name = name;
		this.isDiscrete = isDiscrete;
		this.type = type;
		if (properties != null) propertyMap = new HashMap(properties);
	}
	
	/**
	 * Constructor
	 * @param name			this trait
	 * @param isDiscrete	true if this trait is discrete, false if it is continuous
	 * @param type			the trait type; data, covariate, factor, marker, or exclude
	 */
	public Trait(String name, boolean isDiscrete, String type) {
		this(name, isDiscrete, type, null);
	}
	
	/**
	 * @param trait the trait to be copied
	 * @return a copy of the trait
	 */
	public static Trait getInstance(Trait trait) {
		Trait newTrait = new Trait(trait.name, trait.isDiscrete, trait.type);
		if (trait.propertyMap != null) newTrait.propertyMap = new HashMap<String, Object>(trait.propertyMap);
		return newTrait;
	}
	
	/**
	 * @return	the name by which the trait is identified
	 */
	public String getName() {return name;}

	/**
	 * @return true if this trait is discrete, false if it is continuous
	 */
	public boolean isDiscrete() {return isDiscrete;}
	
	/**
	 * @param levelLabels for a discrete trait, the names of the levels or values this trait can take
	 */
	public void setLevelLabels(String[] levelLabels) {
		setProperty(Trait.PROP_LEVEL_LABELS, levelLabels);
	}
	
	/**
	 * @return for a discrete trait, the names of the levels or values this trait can take 
	 */
	public String[] getLevelLabels() {
        if (propertyMap == null) return null;
		return (String[]) propertyMap.get(PROP_LEVEL_LABELS);
	}
	
	/**
	 * @param level the level number (0,1,...)
	 * @return the name of this level
	 */
	public String getLevelLabel(int level) {
        if (propertyMap == null) return null;
		return getLevelLabels()[level];
	}

	public int getNumberOfLevels() {
        if (propertyMap == null) return 0;
		return getLevelLabels().length;
	}

	/**
	 * @return the number of entries in the property map
	 */
	public int getNumberOfProperties() {
        if (propertyMap == null) return 0;
		return propertyMap.size();
	}

	/**
	 * @return the entry Set for this traits property map
	 */
	public Set<Entry<String, Object>> getProperties() {
        if (propertyMap == null) return null;
		return propertyMap.entrySet();
	}

	/**
	 * @param propertyName the name of a property for this trait
	 * @return the value of this property, null if the property does not exist.
	 */
	public Object getProperty(String propertyName) {
		if (propertyMap == null) return null;
		return propertyMap.get(propertyName);
	}

	/**
	 * @return all of the propery names for this trait
	 */
	public Set<String> getPropertyNames() {
        if (propertyMap == null) return null;
		return propertyMap.keySet();
	}

	/**
	 * Adds a property or sets a new value for an existing property
	 * @param propertyName	a name
	 * @param value			the value for this property
	 */
	public void setProperty(String propertyName, Object value) {
		if (propertyMap == null) propertyMap = new HashMap<String, Object>(); 
		propertyMap.put(propertyName, value);
	}

	/**
	 * Adds a new factor and its value or changes the value of an existing factor. A factor
	 * might be an environment or a rep number.
	 * @param name		this factor's name
	 * @param value		the factor's value for this trait
	 */
	public void addFactor(String name, String value) {
		ArrayList<String> factors = getFactorNames();
		if (factors == null) {
			factors = new ArrayList<String>();
			setProperty(PROP_FACTORS, factors);
		}
		factors.add(name);
		setProperty(name,value);
	}

	/**
	 * @return the number of factors for this trait
	 */
	public int getNumberOfFactors() {
		ArrayList<String> names = getFactorNames();
		if (names == null) return 0;
		return names.size();
	}
	
	/**
	 * @return the names of the factors for this trait. Returns null if there are no factors.
	 */
	public ArrayList<String> getFactorNames() {
		try{return (ArrayList<String>) propertyMap.get(PROP_FACTORS);}
		catch(Exception e){return null;}
	}
	
	/**
	 * @param	factorName the name of a factor
	 * @return	the value of this factor or null if the factor does not exist
	 */
	public String getFactorValue(String factorName) {
		return (String) getProperty(factorName);
	}
	
	/**
	 * @return the type of this trait: data, covariate, factor, marker, or exclude
	 */
	public String getType() {
		return type;
	}

	/**
	 * @param type the type of this trait
	 */
	public void setType(String type) {
		this.type = type;
	}

	public boolean equals(Object obj) {
		if (!(obj instanceof Trait)) return false;
		return toString().equals(obj.toString());
	}

	public int compareTo(Trait trait) {
		return toString().compareTo(trait.toString());
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder(name);
		if (getNumberOfFactors() > 0) {
			ArrayList<String> factorNames = getFactorNames();
			for (String st : factorNames) sb.append(":").append(getFactorValue(st));
		}
		return sb.toString();
	}

	public int hashCode() {
		return toString().hashCode();
	}
	
	public boolean hasLevels() {
		if (propertyMap == null) return false;
		return (propertyMap.get(PROP_LEVEL_LABELS) != null);
	}

}
