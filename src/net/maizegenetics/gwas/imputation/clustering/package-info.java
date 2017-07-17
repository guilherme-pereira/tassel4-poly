/**
 * The clustering package provides classes to cluster genotype sequence for the purpose of identifying common haplotypes. 
 * Missing data creates challenges to clustering that have been addressed by creating clusters of all haplotypes whose pairwise distances are all zero.
 * Because of missing data an individual genotype can belong to more than one haplotype cluster. Resulting clusters can be merged and cluster distances calculated using different methods.
 */

package net.maizegenetics.gwas.imputation.clustering;