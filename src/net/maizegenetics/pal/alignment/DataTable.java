package net.maizegenetics.pal.alignment;

/**
 * Presents data from (at most) one Alignment and one Phenotype in a manner that is 
 * tailored to a set of analyses or display classes. Because DataTable implementations
 * can vary considerably in how they present the underlying data, the main purpose 
 * of the interface is to provide access to the underlying data so that a 
 * new DataTable can be constructed from a different DataTable object. For example
 * a DataTable designed for the Alignment Viewer may be attached to the DataTree. 
 * For a specific analysis, a new DataTable specific to that analysis would be instantiated.
 * @author Peter
 *
 */
public interface DataTable {
	/**
	 * @return	the Alignment which this DataTable contains
	 */
	Alignment getAlignment();
	
	/**
	 * @return	the Phenotype which this DataTable contains
	 */
	Phenotype getPhenotype();
	
	/**
	 * @return	true if this DataTable represents all taxa contained in either 
	 * the Alignment or the Phenotype, false if it represents only those taxa 
	 * contained in both the Alignment and the Phenotype. 
	 */
	boolean isUnion();
}
