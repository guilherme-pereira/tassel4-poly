package net.maizegenetics.pal.alignment;

import java.io.Serializable;
import java.util.List;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.report.Report;
import net.maizegenetics.pal.report.TableReport;

/**
 * A collection of data for a set of Taxa. Each taxon is expected to appear only once in the data. 
 * Data is represented as a two dimensional table of doubles, with taxa as rows and traits as columns. 
 * @author Peter Bradbury
 *
 */
public interface Phenotype extends Serializable, Report, TableReport {
	public static final String MISSING_STRING = "?";
	public static final double MISSING = Double.NaN;
	
    /**
     * @param taxon an integer i, representing the ith row in the data set
     * @param trait an integer j, representing the jth column in the data set
     * @return 		the data stored for the ith taxon, jth trait
     */
    public double getData(int taxon, int trait);

    /**
     * @param taxon
     * @param trait
     * @return		the data value for this taxon and trait
     */
    public double getData(Identifier taxon, Trait trait);

    /**
     * Sets the double value for this trait and taxon.
     * @param taxon	an integer i, representing the ith row in the data set
     * @param trait	an integer j, representing the jth column in the data set
     * @param value	the data to be stored for the ith taxon, jth trait
     */
    public void setData(int taxon, int trait, double value);

    /**
     * Sets the double value for this trait and taxon.
     * @param taxon
     * @param trait
     * @param value
     */
    public void setData(Identifier taxon, Trait trait, double value);

    /**
     * @param trait
     * @return 	the index of the trait, representing the column in which the data for this trait is stored.
     */
    public int whichTrait(Trait trait);

    /**
     * @param taxon
     * @return	the index of the taxon, representing the column in which the data for this trait is stored.
     */
    public int whichTaxon(Identifier taxon);

    /**
     * @param trait	an integer j, representing the jth column in the data set
     * @return		the Trait for this column
     */
    public Trait getTrait(int trait);

    /**
     * @param taxon	an integer i, representing the ith row in the data set
     * @return		the taxon Identifier for this row
     */
    public Identifier getTaxon(int taxon);

    /**
     * @return the number of taxa or rows in this data set
     */
    public int getNumberOfTaxa();

    /**
     * @return the number of traits or columns in this data set
     */
    public int getNumberOfTraits();

    /**
     * @return the traits or columns of this dataset, in order by column number
     */
    public List<Trait> getTraits();

    /**
     * @return	the taxa in this data set in order by row number
     */
    public IdGroup getTaxa();

    /**
     * @return the data stored in this data set as an array of doubles
     */
    public double[][] getData();
    
    //these functions are necessary for backward compatability

    /**
    *
    * @param factor
    * @return the name of this factor
    */
   String getFactorName(int factor);

   /**
    * Use a copy of the factor names to create a new Phenotype to avoid making unintended changes to the original
    * @return
    */
   String[] getFactorNameCopy();


   /**
    *
    * @return
    */
   int getNumberOfFactors();

   
   /**
    * for a factor, sets the active value, which determines whether the factor is to be used in analyses
    * @param factor  the factor number
    * @param active  is the factor to be used in analyses
    */
   void setActiveFactor(int factor, boolean active);
   
   /**
    * is this factor to be used in analyses?
    * @param factor
    * @return  true or false, is factor to be used in analyses
    */
   boolean isFactorActive(int factor);
   
   /**
    * 
    * @param factor this factor
    * @return  the name of the factor
    */
   String getActiveFactorName(int factor);
   
   /**
    * An array of the names of all the active factors
    * @return
    */
   String[] getActiveFactorNames();
   
   /**
    * 
    * @return the number of active factors
    */
   int getNumberOfActiveFactors();
   
   String[] getRowNames();
}
