/*
 * MutableAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.Identifier;

/**
 *
 * @author terry, Gabriel Rodrigues Alves Margarido
 */
public interface MutableAlignment extends Alignment {

    /**
     * Sets base at given taxon and site. newBase should contain both diploid
     * allele values (one in first four bits and second in last four bits). This
     * method may not be supported by every implementation.
     *
     * @param taxon taxon
     * @param site site
     * @param newBase new diploid allele values
     */
    public void setBase(int taxon, int site, byte newBase);

    /**
     * Sets bases starting at specified site for given taxon. Each base should
     * contain both diploid allele values (one in first four bits and second in
     * last four bits). This method may not be supported by every
     * implementation.
     *
     * @param taxon taxon
     * @param startSite starting site
     * @param newBases new diploid allele values
     */
    public void setBaseRange(int taxon, int startSite, byte[] newBases);

    public void addSite(int site);

    public void removeSite(int site);

    /**
     * This clears all data from given site, and when next clean() is performed,
     * it will be removed.
     *
     * @param site site
     */
    public void clearSiteForRemoval(int site);

    /**
     * Adds given identifier (taxon) to the end of taxa list.
     *
     * @param id identifier to add
     */
    public void addTaxon(Identifier id);

    /**
     * Sets the identifier at given taxon index.
     *
     * @param taxon taxon index
     * @param id identifier to set
     */
    public void setTaxonName(int taxon, Identifier id);

    public void removeTaxon(int taxon);

    public void setPositionOfSite(int site, int position);

    public void setLocusOfSite(int site, Locus locus);
    
    /**
     * Sets depth count for each diploid allele at the given taxon and site.
     * 
     * @param taxon taxon
     * @param site site
     * @param values values
     */
    public void setDepthForAlleles(int taxon, int site, short[] values);

    /**
     * Sets possible alleles for a site
     * 
     * @param site
     * @param values 
     */
    public void setCommonAlleles(int site, byte[] values);
    
    /**
     * Clean alignment including sorting sites by position.
     */
    public void clean();

    /**
     * Sets the reference allele for a site
     * 
     * @param site
     * @param diploidAllele
     */
    public void setReferenceAllele(int site, byte diploidAllele);
    
    public void setSNPID(int site, String name);
    
    /**
     * True if changes since last clean().
     *
     * @return true if dirty.
     */
    public boolean isDirty();
}
