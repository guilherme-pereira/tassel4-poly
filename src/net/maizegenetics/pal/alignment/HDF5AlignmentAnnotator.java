/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.alignment;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import java.util.Arrays;

/**
 * Class that provide methods for calculating many of the key computed statics
 * for an HDF5 nucleotide alignment based on TASSEL half bytes.
 *
 * It is designed to be used in multi-thread situations, where one tread may calculate
 * allele frequency while another is calculating allele presence, or bit alignments.
 * 
 * Essentially, this class is called when the HDF5 alignment become dirty to recalculate 
 * all the key statistics.
 * 
 * @author edbuckler
 */
public class HDF5AlignmentAnnotator implements Runnable {
    public static enum AnnoType {ALLELEFreq, SBIT, TBIT};
    private MutableNucleotideAlignmentHDF5 ma;
    private String hdf5File;
    private AnnoType theAnno;

    public HDF5AlignmentAnnotator(MutableNucleotideAlignmentHDF5 ma, String hdf5File, AnnoType theAnno) {
        this.ma=ma;
        this.hdf5File=hdf5File;
        this.theAnno=theAnno;
    }
    
    @Override
    public void run() {
        switch (theAnno) {
            case ALLELEFreq:
                calculateAlleleFreq();
                break;
            default:
                throw new AssertionError();
        }
    }
    
    private void calculateAlleleFreq() {
        IHDF5Writer reader = HDF5Factory.open(hdf5File);
        int sites=ma.getSiteCount();
        int taxa=ma.getSequenceCount();
        int[][] af=new int[6][sites];
        byte[][] afOrder=new byte[6][sites];
        float[] coverage=new float[ma.getSequenceCount()];
        float[] hets=new float[ma.getSequenceCount()];
        for (int taxon = 0; taxon < taxa; taxon++) {
            byte[] genotype=reader.readByteArray(ma.getTaxaGenoPath(taxon));
            int covSum=0;  //coverage of the taxon
            int hetSum=0;
            for (int s = 0; s < sites; s++) {
                byte[] b = AlignmentUtils.getDiploidValues(genotype[s]);
                if(b[0]<6) af[b[0]][s]++;
                if(b[1]<6) af[b[1]][s]++;
                if(AlignmentUtils.isHeterozygous(genotype[s])) hetSum++;
                if(genotype[s]!=Alignment.UNKNOWN_DIPLOID_ALLELE) covSum++;
            }
            coverage[taxon]=(float)covSum/(float)sites;
            hets[taxon]=(float)hetSum/(float)covSum;
        }
//        byte[] mjAlleles=new byte[sites];
//        byte[] mnAlleles=new byte[sites];
        float[] maf=new float[sites];
        float[] paf=new float[sites];        
        int baseMask=0xF;
        for (int s = 0; s < sites; s++) { 
            int sum=0;
            int[] cntAndAllele=new int[6];
            for (byte i = 0; i < 6; i++) {
                cntAndAllele[i]=(af[i][s]<<4)|(5-i);  //size | allele (the 5-i is to get the sort right, so if case of ties A is first)
                sum+=af[i][s];
            }
            Arrays.sort(cntAndAllele);  //ascending quick sort, there are faster ways
            //http://stackoverflow.com/questions/2786899/fastest-sort-of-fixed-length-6-int-array
            for (byte i = 0; i < 6; i++) {
                afOrder[5-i][s]=(cntAndAllele[i]>0xF)?((byte)(5-(baseMask&cntAndAllele[i]))):Alignment.UNKNOWN_ALLELE;
            }
            if(afOrder[1][s]!=Alignment.UNKNOWN_ALLELE) maf[s]=(float)af[afOrder[1][s]][s]/(float)sum;
            paf[s]=(float)sum/(float)(2*taxa);
        } 
        ma.setCalcAlleleFreq(af, afOrder, maf, paf, coverage, hets);
    }
    
    
}
