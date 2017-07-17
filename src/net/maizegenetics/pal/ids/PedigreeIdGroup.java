/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.ids;

/**
 * Index group of identifiers that includes information on relationships between
 * taxa.
 * @author edbuckler
 */
public interface PedigreeIdGroup extends IdGroup {
    /*
     * Returns the parents of a given taxon.  Null if unknown.  Generally the array
     * is size 1 or 2, but if the taxon is derived from a synthetic it could be larger.
     */
    public int[] getParentIndices(int taxon);
    
    /*
     * Returns the children of a given taxon.  Null if unknown.  Generally the array
     * is size 1 or 2, but if the taxon is derived from a synthetic it could be larger.
     */
    public int[] getChildrenIndices(int taxon);
    
}
