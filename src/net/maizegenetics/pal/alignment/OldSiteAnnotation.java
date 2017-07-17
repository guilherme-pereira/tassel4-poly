package net.maizegenetics.pal.alignment;
 
/**
 * Container class for holding information about the site.  This includes information
 * on position, MAF, coverage, 
 * @deprecated 
 * @author edbuckler
 */
@Deprecated
public class OldSiteAnnotation {
    public int position;  //4
    public byte[] myAlleleFreqOrder;  //2-6 bytes, sorted by allele frequency
    public int[] myAlleleCnt;  //2-6*4=24,  sorted by allele frequency of myAlleleFreqOrder
    public float maf; //4
    public float siteCov;  //4
    public String mySNPIDs;  //~14 perhaps don't keep in memory, consider compressed
    //add: byte:reference allele
    //add: byte:number of alleles
    
    public OldSiteAnnotation(int position, byte[] myAlleleFreqOrder, int[] myAlleleCnt, float maf, float siteCov, String mySNPIDs) {
        this.position = position;
        this.myAlleleFreqOrder = myAlleleFreqOrder;
        this.myAlleleCnt = myAlleleCnt;
        this.maf = maf;
        this.siteCov = siteCov;
        this.mySNPIDs = mySNPIDs;
    }
 
    public OldSiteAnnotation(int position) {
        this.position = position;
        this.myAlleleFreqOrder = null;
        this.myAlleleCnt = null;
        this.maf = Float.NaN;
        this.siteCov = Float.NaN;
        this.mySNPIDs = null;
    }
    
    public int[][] getAllelesSortedByFrequency() {
        int result[][] = new int[2][myAlleleCnt.length];
            for (int i = 0; i < myAlleleFreqOrder.length; i++) {
               result[0][i]=myAlleleFreqOrder[i];
               result[1][i]=myAlleleCnt[i];
           }
         return result;
    }
    
    public int getAlleleTotal() {
        int total=0;
        for (int b : myAlleleCnt) {total+=b;}
        return total;
    }
    
    
}
