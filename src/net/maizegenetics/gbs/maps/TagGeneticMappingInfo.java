/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.maps;

/**
 * Class used to to store genetic mapping for hypothesis and output to H5 format
 * @author Fei Lu
 */
public class TagGeneticMappingInfo {
    /**P-value from binomial test*/
    public double p = Double.NEGATIVE_INFINITY;
    /**Chromosome of the most significant site*/
    public int chromosome = Integer.MIN_VALUE;
    /**Position of the most significant site*/
    public int position = Integer.MIN_VALUE;
    /**Number of sites passing P-value threshold*/
    public int sigSiteNum =Integer.MIN_VALUE;
    /**Range of sites passing P-value threshold, for the first significant site to the last*/
    public int sigSiteRange = Integer.MIN_VALUE;
    
    public TagGeneticMappingInfo () {
        
    }
    
    public TagGeneticMappingInfo (double p, int chromosome, int position, int sigSiteNum, int sigSiteRange) {
        this.p = p;
        this.chromosome = chromosome;
        this.position = position;
        this.sigSiteNum = sigSiteNum;
        this.sigSiteRange = sigSiteRange;
    }
    
}
